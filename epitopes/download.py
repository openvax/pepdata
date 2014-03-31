import gzip
from os import (environ, makedirs, remove)
from os.path import (join, exists, isdir, splitext)
from shutil import move, rmtree, copyfileobj
from tempfile import NamedTemporaryFile
import zipfile
import logging
import sqlite3
try:
    from urllib2 import urlretrieve, urlopen
except ImportError:
    from urllib import urlretrieve, urlopen

import appdirs
import pandas as pd
from progressbar import ProgressBar
from Bio import SeqIO

def ensure_dir(path):
    if not exists(path):
        makedirs(path)

def get_data_dir(subdir = None, envkey =  None):
    if subdir is None: subdir = "epitopes"

    dir = appdirs.user_cache_dir(subdir)

    if subdir == "epitopes" and "EPITOPES_data_dir" in environ:
        return environ["EPITOPES_data_dir"]
    else:
        return dir

def fetch_data(filename, download_url, subdir = None):
    data_dir = get_data_dir(subdir)
    ensure_dir(data_dir)
    full_path = join(data_dir, filename)
    if not exists(full_path):
        logging.info("Downloading %s",  download_url)

        base_name, ext = splitext(filename)
        tmp_file = NamedTemporaryFile(
            suffix='.' + ext,
            prefix = base_name,
            delete = False)
        tmp_path = tmp_file.name
        in_stream = urlopen(download_url)
        copyfileobj(in_stream, tmp_file)
        in_stream.close()
        tmp_file.close()

        if download_url.endswith("zip") and not filename.endswith("zip"):
            logging.info("Decompressing zip into %s...", filename)
            with zipfile.ZipFile(tmp_path) as z:
                extract_path = z.extract(filename)
            move(extract_path, full_path)
            remove(tmp_path)
        elif download_url.endswith("gz") and not filename.endswith("gz"):
            logging.info("Decompressing gzip into %s...", filename)
            with gzip.GzipFile(tmp_path) as src:
                contents = src.read()
            remove(tmp_path)
            with open(full_path, 'w') as dst:
                dst.write(contents)
        elif download_url.endswith(("html", "htm")):
            logging.info("Extracting HTML table into CSV %s...", filename)
            df = pd.read_html(tmp_path, header=0, infer_types=False)[0]
            df.to_csv(full_path, sep=',', index=False, encoding='utf-8')
        else:
            move(tmp_path, full_path)

    return full_path

def fetch_fasta_dict(filename, download_url, subdir = None):
   fasta_path = fetch_data(filename, download_url, subdir)
   return SeqIO.index(fasta_path, 'fasta')

def _db_table_exists(db, table_name):
    query = \
        "SELECT name FROM sqlite_master WHERE type='table' AND name='%s'" % \
        table_name
    for _ in db.execute(query):
        return True
    return False

def fetch_fasta_db(db_filename, table_name,
        fasta_filename, download_url, subdir = None):
    data_dir = get_data_dir(subdir)
    db_path = join(data_dir, db_filename)
    db = sqlite3.connect(db_path)
    if _db_table_exists(db, table_name):
        return db
    create = \
        "create table if not exists %s (id TEXT, seq TEXT)" % table_name
    db.execute(create)
    dict = fetch_fasta_dict(fasta_filename, download_url, subdir)
    rows = [
        (idx, str(record.seq))
        for (idx, record)
        in dict.iteritems()
    ]
    db.executemany("insert into %s values (?, ?)" % table_name, rows)
    db.commit()
    return db


def fetch_and_transform_data(
        transformed_filename,
        transformer,
        loader,
        source_filename,
        source_url,
        subdir = None):
    data_dir = get_data_dir(subdir)
    ensure_dir(data_dir)
    transformed_path = join(data_dir, transformed_filename)
    if not exists(transformed_path):
        source_path = fetch_data(source_filename, source_url, subdir)
        result = transformer(source_path, transformed_path)
    else:
        result = loader(transformed_path)
    assert exists(transformed_path)
    return result


def clear_cache(subdir = None):
    data_dir = get_data_dir(subdir)
    rmtree(data_dir)
