# encoding: utf-8
#
# @Author: Tom Dwelly
# @Date: Late 2019
# @Filename: catalogdb_wrapper.py
# @License: BSD 3-Clause
# @Copyright: Tom Dwelly
# Content:
#    wrapper to aid with remote connections to the sdss5db

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import socket
import nclib

from sdssdb.connection import PeeweeDatabaseConnection

from sdssdb.peewee.sdss5db import database
from sdssdb.peewee.sdss5db import catalogdb

import numpy as np
import pandas as pd

from .print_func import *

def setup_db_connection():
    '''
     Sets up a remote connections to the sdss5db, using a SSH tunnel if it exists
    '''

    #sdss5db_host = 'sdssadmin.wasatch.peaks'
    sdss5db_host = 'operations-test.sdss.utah.edu'

    # test if we are already connected
    if database.connected is True:
        return True

    database.autorollback = True
    database.set_profile('sdssadmin')

    hostname = socket.gethostname()
    # determine which machine we are working on
    if hostname == 'eboss':
        database.connect_from_parameters(user='sdss', host=sdss5db_host, port=5432)

    else:
        envname = 'LOCAL_PORT_UTAH_SDSS5DB'
        sdss5db_port = os.getenv(envname)
        if sdss5db_port is None :
            msg = f'Error, you must set the following environment variable: {envname}'
            raise Exception(msg)
        else:
            sdss5db_port = int(sdss5db_port)

        # test for an existing ssh tunnel
        try:
            nc = nclib.Netcat(('localhost', sdss5db_port))
            print_comment(f'Using existing ssh tunnel to DB server via localhost:{sdss5db_port}')
            del nc
            database.connect_from_parameters(user='sdss', host='localhost', port=sdss5db_port)
        except:
            print_error(f'You first need to set up an SSH tunnel to sdssdb on port: {sdss5db_port}')
            print_error(f'Hint: ssh -l your_utah_username -L {sdss5db_port}:{sdss5db_host}:5432 eboss.sdss.org cat -')

    return database.connected



def retrieve_targets_from_catalogdb(tablename=None,
                                    constraints=[],
                                    fmt='pandas',
                                    columns=['ra',
                                             'dec'],
                                    ):

    '''Query objects in a single table, apply the given list
    of constraints, retrieving the given subset of columns,
    and returning the result in the given data format.
    '''

    assert tablename is not None

    # try to connect to the database
    if setup_db_connection() is False :
        raise Exception('Unable to connect to database')

    modelname =  f'catalogdb.{tablename}'
    try :
        db = eval(modelname)
    except:
        raise Exception(f'Problem initialising database model: {modelname}')


    ## turn the list of columns as strings into a list of model fields
    collist = []
    for col in columns:
        try :
            collist.append(eval(f'db.{col}'))
        except:
            raise Exception(f'Specified field ({col}) not found in model {modelname}')

    targets_query = db.select(*collist)

    # apply the constraints in series
    if constraints:
        for c in constraints:
            targets_query = targets_query.where(eval(c))

    print_comment(f'This selection will return {targets_query.count()} objects')
    print_comment(f'Retrieving objects from DB into Pandas DataFrame')

    try:
        df = pd.DataFrame(list(targets_query.dicts()))
    except:
        raise Exception(f'Database query failed: {targets_query}')

    if ( fmt is None or fmt == '' ):
        return df
    elif fmt.lower() in ['pandas', 'pd', 'df'] :
        return df
    elif fmt.lower() in ['numpy', 'np', 'recarray'] :
        print_comment(f'Converting result to numpy recarray')
        return df.to_records()
    else:
        return None
