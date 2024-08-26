
import io
import os
import re
import json
import base64


from jinja2 import PackageLoader,Environment,FileSystemLoader
from utils import lzstring

def compress_json(data):
    """Take a Python data object. Convert to JSON and compress using lzstring"""
    json_string = json.dumps(data).encode("utf-8", "ignore").decode("utf-8")
    json_string = sanitise_json(json_string)
    x = lzstring.LZString()
    return x.compressToBase64(json_string)

def sanitise_json(json_string):
    """
    The Python json module uses a bunch of values which are valid JavaScript
    but invalid JSON. These crash the browser when parsing the JSON.
    Nothing in the MultiQC front-end uses these values, so instead we just
    do a find-and-replace for them and switch them with `null`, which works fine.

    Side effect: Any string values that include the word "Infinity"
    (case-sensitive) will have it switched for "null". Hopefully that doesn't happen
    a lot, otherwise we'll have to do this in a more complicated manner.
    """
    json_string = re.sub(r"\bNaN\b", "null", json_string)
    json_string = re.sub(r"\b-?Infinity\b", "null", json_string)
    return json_string


# from .utils import config
report = dict()

config = dict()
config['title'] = 'bsDoctor'


# data = [1, 'foo']
data2 = [1,2,['a', 3]]

# data = [{'x':x, 'y':y} for x,y in zip(range(10),range(10))]
# data = [list(range(10)), list(range(5))]

data = {
    'x': list(range(10)),
    'y': list(range(10))
}

## compressed json data
# report['plot_compressed_json'] = compress_json(data)
report['plot_compressed_json'] = data

# Function to include file contents in Jinja template
tmp_dir = 'report/'
def include_file(name, fdir=tmp_dir, b64=False):
    try:
        if fdir is None:
            fdir = ""
        if b64:
            with io.open(os.path.join(fdir, name), "rb") as f:
                return base64.b64encode(f.read()).decode("utf-8")
        else:
            with io.open(os.path.join(fdir, name), "r", encoding="utf-8") as f:
                return f.read()
    except (OSError, IOError) as e:
        # logger.error(f"Could not include file '{name}': {e}")
        pass

env = Environment(loader=FileSystemLoader('report/'))
env.globals["include_file"] = include_file
template = env.get_template('base.html')    
temp_out = template.render(data=data, report=report, config=config)   

with open('report/out.html', 'w', encoding='utf-8') as f:
    f.writelines(temp_out)
    f.close()
