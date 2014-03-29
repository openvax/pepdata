import gzip
from os import (environ, makedirs, remove)
from os.path import (join, exists, isdir, splitext)
from shutil import move, rmtree, copyfileobj
from tempfile import NamedTemporaryFile
import zipfile
try:
    from urllib2 import urlretrieve, urlopen
except ImportError:
    from urllib import urlretrieve, urlopen

import appdirs
import pandas as pd
from progressbar import ProgressBar

def ensure_dir(path):
    if not exists(path):
        makedirs(path)
DATA_DIR = environ.get("EPITOPES_DATA_DIR", appdirs.user_cache_dir("epitopes"))

def fetch_data(filename, download_url):
    ensure_dir(DATA_DIR)
    full_path = join(DATA_DIR, filename)
    if not exists(full_path):
        print "Downloading %s" %  download_url

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

        if download_url.endswith("zip"):
            print "Decompressing..."
            with zipfile.ZipFile(tmp_path) as z:
                extract_path = z.extract(filename)
            move(extract_path, full_path)
            remove(tmp_path)
        elif download_url.endswith("gz"):
            print "Decompressing..."
            with gzip.GzipFile(tmp_path) as src:
                contents = src.read()
            remove(tmp_path)
            with open(full_path, 'w') as dst:
                dst.write(contents)
        elif download_url.endswith(("html", "htm")):
            df = pd.read_html(tmp_path, header=0, infer_types=False)[0]
            df.to_csv(full_path, sep=',', index=False, encoding='utf-8')
        else:
            move(tmp_path, full_path)

    return full_path

def fetch_and_transform_data(
        transformed_filename,
        transformer,
        loader,
        source_filename,
        source_url):
    ensure_dir(DATA_DIR)
    transformed_path = join(DATA_DIR, transformed_filename)
    if not exists(transformed_path):
        source_path = fetch_data(source_filename, source_url)
        result = transformer(source_path, transformed_path)
    else:
        result = loader(transformed_path)
    assert exists(transformed_path)
    return result

def clear_cache():
    rmtree(DATA_DIR)
