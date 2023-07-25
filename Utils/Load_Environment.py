#!/usr/bin/env python3
import sys
import os
import pwd

from Scope.Classes_Input import *
from Scope.Classes_Environment import *
from Scope.Classes_Queue import *

def load_environment(env_path: str, debug: int=0):
    assert type(env_path) == str
    if os.path.isfile(env_path): 
        try:
            env = load_binary(env_path)
        except Exception as exc:
            print(exc)
    return env

def configure_environment(set_user_queues: bool=True, check_usage: bool=True):
    env=environment()
    env.get_mqueues()
    if set_user_queues: 
        env.user_queues_preferences()
        if check_usage: 
            for q in env.queues:
                q.check_usage
    return env

#        if env_path is None: 
#            savefile_path = '/home/'+env.user+'/'
#            savefile_name = 'scope_environment.npy'
#            if os.path.isdir(savefile_path): env.save(savefile_path+savefile_name)
#    
#        elif type(env_path) == str: 
#            env.save(env_path)
#        else: print("LOAD_ENVIRONMENT: I could not understand the savefile_path:", savefile_path)
#    
#    print("ENVIRONMENT LOADED and saved here:", savefile_path) 
#    print(env)
#
#    return env
