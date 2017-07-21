"""
This module contains classes necessary to introduce new data types to the SQLAlchemy
ORM as well as other helper functions.
"""

from __future__ import absolute_import

from sqlalchemy import Table, Column, func, event
from sqlalchemy.schema import DDL
from sqlalchemy.dialects.postgresql import BYTEA
import sqlalchemy.types as types

try:
    import rdkit
except ImportError:
    print "Failed to import RDKit module. Tables and operations dependent on it will have errors."

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
        return "ptree"

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

class RDMol(types.TypeDecorator):
    """Custom type for RDKit Molecule type""" 
    impl = BYTEA

    def get_col_spec(self):
        return "rdkit.mol"

    def process_bind_param(self, value, dialect):
        if value is not None:
            value = value.ToBinary()
        return value
        
    def bind_expression(self, bind_value):
        return func.rdkit.mol_from_pkl(bind_value)

    def process_result_value(self, value, dialect):
        return rdkit.Mol(str(value)) if value is not None else value
        
    def column_expression(self, col):
        return func.rdkit.mol_to_pkl(col, type_=self)


class RDBinaryFingerprint(types.UserDefinedType):
    """Custom type for RDKit Fingerprint type""" 
    def __init__(self):
        pass

    def get_col_spec(self):
        return "rdkit.bfp"

    def bind_processor(self, dialect):
        def process(value):
            return value
        return process

    def result_processor(self, dialect, coltype):
        # \x090000000900888000000707000000807606
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
    event.listen(table, "after_create", CreateComment(table, comments['table']))

    # add comments to the table columns as well
    for column in table.columns:
        try:
            comment = comments['columns'].pop(column.name)
        except KeyError: pass

        event.listen(table, "after_create", CreateComment(column, comment))
        
    if len(comments['columns']):
        print "WARNING: The following provided column comments for table %s were not applied:\n%s" % (table.name, comments['columns'])


def create_partition_insert_trigger(table, interval=25000, column='biomolecule_id', part_label='biomol'):

    """
    Creates a trigger on the given (partitioned) table to insert entries into their correct partition based on biomolecule_id
    """

    create_func = """
        CREATE OR REPLACE FUNCTION %(schema)s.%(table)s_insert_trigger()
        RETURNS TRIGGER AS $$
        DECLARE
            intvl integer;
            part_bound integer;
        BEGIN
            intvl = TG_ARGV[0]::integer;
            IF NEW.{col_id} %% intvl = 0 THEN
                part_bound := NEW.{col_id};
            ELSE
                part_bound := intvl * ((NEW.{col_id} / intvl) + 1);
            END IF;
            EXECUTE format('INSERT INTO %(fullname)s_{part_label}_le_%%s VALUES ($1.*);', part_bound) USING NEW;
            RETURN NULL;
        END;
        $$
        LANGUAGE plpgsql;
        """.format(col_id=column, part_label=part_label)

    create_trigger = """
        CREATE TRIGGER insert_%(table)s_trigger
        BEFORE INSERT ON %(fullname)s
        FOR EACH ROW EXECUTE PROCEDURE %(schema)s.%(table)s_insert_trigger(%(intvl)s);
        """

    event.listen(table, 'after_create', DDL(create_func))
    event.listen(table, 'after_create', DDL(create_trigger, context={'intvl': interval}))
