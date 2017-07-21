require "rubygems"
require "logger"
require "date"
require "fileutils"
require "pathname"
require "net/smtp"

include FileUtils

$config = {
  :mmcif_mirror_dir         => Pathname.new("/tlbnas/mirror/pdb/data/structures/all/mmCIF"),
  :schema_map_file          => Pathname.new("/usr/local/db-loader-v4.2.0/db-loader/test/schema_map_pdbx_na.cif"),
  :db_loader_bin            => Pathname.new("/usr/local/db-loader-v4.2.0/bin/db-loader"),
  :db_host                  => "localhost",
  :db_name                  => "MMCIF_DEV",
  :db_dbms                  => "mysql",
  :db_user                  => "bernardo",
  :db_pass                  => "",
  :db_manager               => "bernardo",
  :my_email                 => "bernardo@cryst.bioc.cam.ac.uk",
  :gloria_email             => "dohaeris@cryst.bioc.cam.ac.uk",
  :schema_load_sql_file     => "DB_LOADER_SCHEMA.sql",
  #:schema_load_mod_sql_file => "DB_LOADER_SCHEMA_MOD.sql",
  :schema_drop_sql_file     => "DB_LOADER_SCHEMA_DROP.sql",
  :data_load_sql_file       => "DB_LOADER_LOAD.sql",
  #:field_delimiter          => "'@\\t'",
  :field_delimiter          => "'\\t'",
  #:row_delimiter            => "'#\\n'",
  :row_delimiter            => "'\\n'",
  :temp_dir                 => Pathname.new("/spunky/temp/mmCIF"),
  :log_file                 => "mmcif_rake.log",
  :nthreads                 => 4
}

#$logger_formatter = Logger::Formatter.new

#class Logger
#  def format_message(severity, timestamp, progname, msg)
#    $logger_formatter.call(severity, timestamp, progname, msg)
#  end
#end

class MultiIO
  def initialize(*targets)
    @targets = targets
  end

  def write(*args)
    @targets.each {|t| t.write(*args) }
  end

  def close
    @targets.each(&:close)
  end
end

#$log_file     = STDOUT
$logger       = Logger.new MultiIO.new(STDOUT, File.open($config[:log_file], 'a'))
$logger.level = Logger::DEBUG

def send_email(from       = $config[:my_email],
               from_alias = 'Alicia Higueruelo',
               to         = $config[:gloria_email],
               to_alias   = 'Gloria',
               subject    = "[Gloria] MMCIF on spunky updated",
               message    = "MMCIF database has been updated.")
  msg = <<END_OF_MESSAGE
From: #{from_alias} <#{from}>
To: #{to_alias} <#{to}>
Subject: #{subject}

#{message}
END_OF_MESSAGE

  Net::SMTP.start('localhost') do |smtp|
    smtp.send_message msg, from, to
  end
end

$logger.debug "Initializing..."
# Tasks

desc "Simply just build MMCIF database!"
task :default => [
  #"check:week",
  #"prepare:temp_dir",
  :req_passwd,
  "prepare:files",
  "create:list",
  "create:revised_schema_mapping",
  "create:updated_schema_mapping",
  "create:schema",
  "create:dumps",
  "drop:tables",
  "create:tables",
  "modify:load_sql",
  "import:dumps",
  #"cleanup:tmpdir",
  "send:email"
]

desc "Request SQL password"
task :req_passwd do
  `stty -echo`
  print "Input SQL Password for user '#{$config[:db_user]}': "
  $config[:db_pass] = $stdin.gets.chomp
  `stty echo`
  puts ""
end

namespace :check do

  desc "Check if this week is for MMCIF"
  task :week do
    if Date.today.cweek % 2 == 0
      $logger.debug "This week is not mine."
      puts "This week is not mine."
      exit
    else
      $logger.debug "This week is mine"
      puts "This week is mine."
    end
  end
end


namespace :prepare do

  desc "Prepare a scratch directory"
  task :temp_dir do

    dir = $config[:temp_dir]

    if File.exists? dir
      rm_rf dir, :secure => true
      $logger.info "Removing #{dir}: done"
    end

    mkdir_p dir

    $logger.info "Creating #{dir}: done"
  end


  desc "Uncompress and copy mmCIF files to working directory"
  task :files do

    zipped_files = Dir[$config[:mmcif_mirror_dir].join("*.gz").to_s]
    pidhash = Hash.new

    zipped_files.each_with_index do |zipped_file, i|
      unzipped_file = $config[:temp_dir].join(File.basename(zipped_file, ".gz")).to_s
      if File.size? unzipped_file
        $logger.debug "File #{unzipped_file} already exists. Skipping"
        next
      end

      $logger.debug "Processing #{zipped_file} > #{unzipped_file}"

      pidhash[spawn("gzip -cd #{zipped_file}", :out=>unzipped_file)] = [ zipped_file, unzipped_file ]
      #system "gzip -cd #{zipped_file} > #{unzipped_file}"

      while pidhash.size >= $config[:nthreads] do
        finished, status = Process.wait2
        $logger.debug "Uncompressed '#{pidhash[finished][0]}' to '#{pidhash[finished][1]}': done (#{i+1}/#{zipped_files.size})"
        #$logger.debug "Uncompressing '#{zipped_file}' to '#{unzipped_file}': done (#{i+1}/#{zipped_files.size})"
        pidhash.delete finished
      end
    end
    Process.waitall
    $logger.info "Uncompressing #{zipped_files.size} PDB mmCIF files to #{$config[:temp_dir]}: done"
  end
end


namespace :create do
  desc "Create a LIST file"
  task :list do
    File.open($config[:temp_dir].join("LIST"), 'w') do |file|
      file.puts Dir[$config[:temp_dir].join("*.cif").to_s].map { |f| File.basename(f) }
    end
    $logger.info "Creating a LIST file: done"
  end

  desc "Create a revised schema mapping file"
  task :revised_schema_mapping do
    $logger.debug "Staring creation of revised schema mapping file: #{$config[:temp_dir]}/revised_schema_mapping.cif"
    Dir.chdir $config[:temp_dir] do
      success = system "#{$config[:db_loader_bin]} " +
                       "-map #{$config[:schema_map_file]} " +
                       "-list LIST " +
                       "-revise revised_schema_mapping.cif"
      if success
        $logger.info "Creating a revised schema mapping file: done"
      else
        $logger.error "Error creating revised schema mapping file!"
        break
      end
    end
  end

  desc "Create an updated schema mapping file"
  task :updated_schema_mapping do
    $logger.debug "Staring creation of updated schema mapping file: #{$config[:temp_dir]}/updated_schema_mapping.cif"
    Dir.chdir $config[:temp_dir] do
      success = system "#{$config[:db_loader_bin]} " +
                       "-map #{$config[:schema_map_file]} " +
                       "-update updated_schema_mapping.cif " +
                       "-revise revised_schema_mapping.cif"
      if success
        $logger.info "Creating an updated schema mapping file: done"
      else
        $logger.error "Error creating updated schema mapping file!"
        break
      end
    end
  end

  desc "Create MMCIF RDB schema"
  task :schema do
    Dir.chdir $config[:temp_dir] do
      success = system "#{$config[:db_loader_bin]} " +
                       #"-map #{$config[:schema_map_file]} " +
                       "-map updated_schema_mapping.cif " +
                       "-server #{$config[:db_dbms]} " +
                       "-db #{$config[:db_name]} " +
                       "-schema"
      if success
        $logger.info "Creating MMCIF schema: done"
      else
        $logger.error "Error creating MMCIF schema!"
        break
      end
    end
  end

  desc "Create MMCIF dump files"
  task :dumps do
    Dir.chdir $config[:temp_dir] do
      success = system "#{$config[:db_loader_bin]} " +
                       #"-map #{$config[:schema_map_file]} " +
                       "-map updated_schema_mapping.cif " +
                       "-server #{$config[:db_dbms]} " +
                       "-db #{$config[:db_name]} " +
                       "-ft #{$config[:field_delimiter]} " +
                       "-rt #{$config[:row_delimiter]} " +
                       "-bcp " +
                       "-list LIST " +
                       "-revise revised_schema_mapping.cif"
      if success
        $logger.info "Creating MMCIF dump files: done"
      else
        $logger.error "Error creating MMCIF dump files!"
        break
      end
    end
  end

  desc "Create MMCIF tables"
  task :tables => [:req_passwd] do
    success = system "mysql -f " +
                     "-h #{$config[:db_host]} " +
                     "-u #{$config[:db_user]} " +
                     "-p#{$config[:db_pass]} " +
                     "< #{$config[:temp_dir].join($config[:schema_load_sql_file])}"
    if success
      $logger.info "Creating MMCIF tables: done"
    else
      $logger.error "Error creating MMCIF tables!"
      break
    end
  end

end

namespace :drop do
  desc "Drop MMCIF tables"
  task :tables => [:req_passwd] do
    success = system "mysql -f " +
                     "-h #{$config[:db_host]} " +
                     "-u #{$config[:db_user]} " +
                     "-p#{$config[:db_pass]} " +
                     "< #{$config[:temp_dir].join($config[:schema_drop_sql_file])}"
    if success
      $logger.info "Dropping MMCIF tables: done"
    else
      $logger.error "Error dropping MMCIF tables!"
      break
    end
  end
end


namespace :modify do
  desc "Modify DB_LOADER_LOAD.sql (to prevent the loading of the massive atom_site table)"
  task :load_sql do
    sql_file      = $config[:temp_dir].join($config[:data_load_sql_file])
    FileUtils.copy_file sql_file, sql_file.sub('.sql', '.sql_orig'), :verbose => true
    load_sql      = IO.readlines(sql_file)
    #atom_site_sql = load_sql.slice!(3..6)
    #load_sql << atom_site_sql

    File.open(sql_file, 'w') do |f|
      skip_flag = false
      load_sql.each do |l|
        $logger.debug l
        if l.strip.end_with? '"atom_site.bcp"' or l.strip.end_with? '"atom_site_anisotrop.bcp"'
          skip_flag = true
        end
        
        f.print "-- " if skip_flag
        f.puts l          
          
        if l.strip.end_with? ';'
          skip_flag = false
        end
      end
      #f.puts load_sql.join
    end
    $logger.info "Finished modifying #{$config[:data_load_sql_file]}"
  end
end


namespace :import do
  desc "Import MMCIF dump files to database"
  task :dumps => [:req_passwd] do
    Dir.chdir $config[:temp_dir] do
      my_cmdline = "mysql --local-infile -v -f " +
                   "-h #{$config[:db_host]} " +
                   "-u #{$config[:db_user]} " +
                   "-p#{$config[:db_pass]} " +
                   "< #{$config[:data_load_sql_file]}"
      #$logger.debug my_cmdline
      success = system my_cmdline
      
      if success
        $logger.info "Importing mmCIF dump files to #{$config[:db_host]}.#{$config[:db_name]}: done"
      else
        $logger.error "Error importing mmCIF dump files to #{$config[:db_host]}.#{$config[:db_name]}!"
        break
      end
    end
  end
end


namespace :cleanup do
  desc "Clean up temporary mmCIF files"
  task :tmpdir do
    rm_rf($config[:temp_dir], :secure => true) if Dir.exists? $config[:temp_dir]
  end
end


namespace :send do

  desc "Send log to Gloria team"
  task :email do
    send_email
  end
end
