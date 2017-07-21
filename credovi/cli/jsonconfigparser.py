import os
import json
from UserDict import UserDict

from cement.core import config

class JSONConfigParserHandler(config.CementConfigHandler):
    class Meta:
        interface = config.IConfig
        label = 'JSONConfigParser'

    options = {}

    def __getitem__(self, item):
        """
        """
        return self.options[item]

    def add_section(self, section):
        """
        Add a new section if it doesn't exist.
        Returns: None
        """
        if section not in self.options: self.options[section] = {}

    def get(self, section, key):
        """
        Return a configuration value based on [section][key]. The return value type
        is unknown.

        Required Arguments:

        section
        The [section] of the configuration to pull key value from.
        key
        The configuration key to get the value from.
        Returns: unknown
        """
        if section in self.options: return self.options[section].get(key)
        else: return None

    def get_section_dict(self, section):
        """
        Return a dict of configuration parameters for [section].

        Required Arguments:

        section
            The config [section] to generate a dict from (using that sections keys).

        Return: dict
        """
        if section in self.options: return self.options[section]
        else: return {}

    def get_sections(self):
        """
        Return a list of configuration sections. These are designated by a [block] label in a config file.

        Returns: list
        """
        return self.options.keys()

    def has_section(self, section):
        """
        Returns whether or not the section exists.

        Returns: bool
        """
        return section in self.options

    def has_key(self, section, key):
        """
        """
        try: return self.options[section].has_key(key)
        except KeyError: return False

    def keys(self, section):
        """
        Return a list of configuration keys from section.

        Required Arguments:

        section
        The config [section] to pull keys from.
        Returns: list
        """
        return self.options[section].keys()

    def merge(self, dict_obj, override=True):
        """
        Merges a dict object into the configuration.

        Required Arguments:

        dict_obj
            The dict to merge into the config
        Optional Arguments:

        override
        Whether to override existing values. Default: True
        """
        for section in dict_obj:
            if isinstance(dict_obj[section], dict):
                if not section in self.get_sections():
                    self.add_section(section)

                for key in dict_obj[section]:
                    if override:
                        self.set(section, key, dict_obj[section][key])
                    else:
                        # only set it if the key doesn't exist
                        if not self.has_key(section, key):
                            self.set(section, key, dict_obj[section][key])

    def parse_file(self, file_path):
        """
        Parse config file settings from file_path. Returns True if the file existed,
        and was parsed successfully. Returns False otherwise.

        Required Arguments:

        file_path
            The path to the config file to parse.

        Returns: boolean
        """
        SUCCESS = False

        if os.path.isfile(file_path):
            try:
                self.options.update(json.loads(open(file_path).read()))
            except StandardError:
                print "Failed to parse JSON: %s" % file_path
                raise

    def set(self, section, key, value):
        """
        Set a configuration value based at [section][key].

        Required Arguments:

        section
            The [section] of the configuration to pull key value from.
        key
            The configuration key to set the value at.
        value
            The value to set.
        """
        self.options[section][key] = value

    def setup(self, defaults):
        """
        The setup function is called during application initialization and must
        'setup' the handler object making it ready for the framework or the application
        to make further calls to it.

        Required Arguments:

        defaults
            The application default config dictionary. This is not a config object,
            but rather a dictionary which should be obvious because the config handler
            implementation is what provides the application config object.

        Returns: n/a
        """
        self.merge(defaults)
