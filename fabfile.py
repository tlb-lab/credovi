from fabric.api import *

env.hosts = ['bahamut.bioc.cam.ac.uk:289']
env.user  = 'adrian'

def remote_info():
    run('uname -a')

def local_info():
    local('uname -a')