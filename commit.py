
import os, re
import subprocess
from contextlib import contextmanager

@contextmanager
def change_directory(directory):
    current_directory = os.getcwd()
    os.chdir(directory)
    try:
        yield
    finally:
        os.chdir(current_directory)

def commit(directory):
    with change_directory(directory):
        subprocess.run(["git", "add", "."])
        subprocess.run(["git", "commit", "-m", "Small corrections"])
        subprocess.run(["git", "push", "origin", "main"])


directories = set()
    
# Iterate over all directories
for directory in os.listdir("."):
    if os.path.isdir(directory):
        for dirname in os.walk(directory):
            directories.add( re.sub( "\(|\'","", str(str(dirname).split(",")[0].split("/")[0]) ) )

for d in directories:
    commit ( d )
