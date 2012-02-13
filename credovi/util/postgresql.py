from subprocess import Popen
from collections import Iterable
from multiprocessing import Pool
from tempfile import NamedTemporaryFile

def _copy(filename, table, server, port, db, user):
    """
    """
    statement = "COPY {} FROM '{}' DELIMITER E'\t';".format(table, filename)
    
    cmd = Popen(['psql', '-h', server, '-p', str(port), '-U', user, '-d', db,
                 '-c', statement])
    cmd.wait()

def _copy_from_file(iterable, table, server, port, db, user, parallel=False):
    """
    This method take a list of filenames as well as connection options. It will
    write a temporary script containing all the necessary COPY commands to load
    the files into the specified table. The temporary script is then executed
    with one psql command. Using this function is the best way to load files in
    a VirtualBox environment because it avoids the crappy networking.
    
    Remember that COPY FROM only looks for files on the server!
    """
    # create a temporary file to which all COPY commands will be written
    flike = NamedTemporaryFile()            
    
    # iterate through all the files in our 
    for filename in iterable:
        
        # create a proper COPY statement for this file
        statement = "COPY {} FROM '{}' DELIMITER E'\t';".format(table, filename)
        
        # and write it to the temp file
        flike.file.write(statement)

    # important: flush first before using psql to execute the script
    flike.file.flush()
    
    # send an SSH command to run GNU parallel
    if parallel:

        # not implemented yet - SSH requires a password        
        pass
    
    else:
        
        # create psql command through Popen
        cmd = Popen(['psql', '-h', server, '-p', str(port), '-U', user, '-d', db,
                     '-f', flike.name])
    
        # run!
        cmd.wait()
    
    # delete the temporary file
    flike.close()

def copy(str_or_iter, table, server='bahamut.bioc.cam.ac.uk', port=5432, db='cryst',
         user='adrian', **kwargs):
    """
    """
    # use multiple processes to copy files into the database;
    # does not work with VirtualBox!
    if kwargs.get('processes') > 1:
        
        # does not make any sense to use multiprocessing with only one file
        if not isinstance(str_or_iter, Iterable):
            raise RuntimeError("cannot use multiprocessing with a single file.")
        
        # load them asynchronously through multiprocessing
        pool = Pool(processes=kwargs.get('processes'))
            
        # iterate through all the filenames
        for filename in str_or_iter:            

            # tuple contains the arguments for the pgload function
            pool.apply_async(_copy, (filename, table, server, port, db, user))
        
        # execute
        pool.close()
        pool.join()
    
    # do not use multiprocessing and write all COPY commands to a single file
    # that will be executed with a single psql command afterwards;
    # this has the advantage of causing less stress on the network
    else:
        if isinstance(str_or_iter, Iterable):
            
            #
            _copy_from_file(str_or_iter, table, server, port, db, user)
        
        # run a simple command for this single file
        else: _copy(filename, table, server, port, db, user)
       