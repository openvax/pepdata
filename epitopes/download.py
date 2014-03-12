from os import (environ, makedirs, remove)
from os.path import (join, exists, isdir)
from shutil import move
import gzip
import zipfile

try:
    from urllib2 import urlretrieve
except ImportError:
    from urllib import urlretrieve

import appdirs
import pandas as pd

DATA_DIR = environ.get("EPITOPES_DATA", appdirs.user_cache_dir("epitopes"))

def fetch_data(filename, download_url):
    if not exists(DATA_DIR):
        makedirs(DATA_DIR)
    full_path = join(DATA_DIR, filename)
    if not exists(full_path):
        print "Downloading %s" % download_url
        if download_url.endswith(("csv", "hdf", "txt")):
            urlretrieve(download_url, full_path)
        elif download_url.endswith("zip"):
            tmp_path, _ = urlretrieve(download_url)
            print "Decompressing..."
            with zipfile.ZipFile(tmp_path) as z:
                extract_path = z.extract(filename)
            move(extract_path, full_path)
            remove(tmp_path)
        elif download_url.endswith("gz"):
            tmp_path, _ = urlretrieve(download_url)
            print "Decompressing..."
            with gzip.GzipFile(tmp_path) as src:
                contents = src.read()
            remove(tmp_path)
            with open(full_path, 'w') as dst:
                dst.write(contents)
        elif download_url.endswith(("html", "htm")):
            df = pd.read_html(download_url, header=0, infer_types=False)[0]
            df.to_csv(full_path, sep=',', index=False, encoding='utf-8')
    return full_path
