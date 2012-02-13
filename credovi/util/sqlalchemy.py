"""
This module contains classes necessary to introduce new data types to the SQLAlchemy
ORM as well as other helper functions.
"""

from __future__ import absolute_import

from sqlalchemy import Table, Column
from sqlalchemy.schema import DDL
from sqlalchemy.event import listen
import sqlalchemy.types as types

class ArrayXi(types.UserDefinedType):
    """
    Custom data type to create a vector composite type.
    """
    def __init__(self):
        pass

    def get_col_spec(self):
        return "arrayxi"

    def bind_processor(self, dialect):
        def process(value):
            return value
        return process

    def result_processor(self, dialect, coltype):
        def process(value):
            return value
        return process

class Cube(types.UserDefinedType):
    """
    Data type representing multidimensional cubes.
    """
    def __init__(self):
        pass

    def get_col_spec(self):
        return "CUBE"

    def bind_processor(self, dialect):
        def process(value):
            return value
        return process

    def result_processor(self, dialect, coltype):
        def process(value):
            return value
        return process
    
class PTree(types.UserDefinedType):
    """
    Custom data type to create a vector composite type.
    """
    def __init__(self):
        pass

    def get_col_spec(self):
        return "ltree"

    def bind_processor(self, dialect):
        def process(value):
            return value
        return process

    def result_processor(self, dialect, coltype):
        def process(value):
            return value
        return process

class Vector3D(types.UserDefinedType):
    """
    Custom data type to create a vector composite type.
    """
    def __init__(self):
        pass

    def get_col_spec(self):
        return "VECTOR3D"

    def bind_processor(self, dialect):
        def process(value):
            return map(float, values[1:-1].split(','))
        return process

    def result_processor(self, dialect, coltype):
        def process(value):
            return value
        return process

def CreateComment(element, comment):
    """
    Returns a DDL to comment on a column. This function can only be used in the
    context of an SQLAlchemy listen event or any other function that supports the
    SQLAlchemy %(fullname)s string formatting.
    """
    if isinstance(element, Table):
        statement = "COMMENT ON TABLE %(fullname)s IS '{0}'".format(comment)
    
    elif isinstance(element, Column):
        statement = "COMMENT ON COLUMN %(fullname)s.{0} IS '{1}'".format(element.name, comment)
        
    return DDL(statement)

def comment_on_table_elements(table, comments):
    """
    This function takes a table and a dictionary of comments in the form
    {"table": "table comment", "columns": {"column name": "column comment"}} and
    adds these comments to the table as well as to all its columns through the
    SQLAlchemy event.listen() framework. 
    """
    # add a comment to describe the table
    listen(table, "after_create", CreateComment(table, comments['table']))
    
    # add comments to the table columns as well
    for column in table.columns:
        try:
            comment = comments['columns'][column.name]
        except KeyError: pass
        
        listen(table, "after_create", CreateComment(column, comment))