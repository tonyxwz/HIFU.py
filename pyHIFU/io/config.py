"""
classes and methods to handle json config files
"""
import json


def readjson(json_str=None, json_path=None):
    if json_path is not None:
        with open(json_path, "rb") as f:
            j = json.load(f)
    elif json_str is not None:
        j = json.loads(json_str)

    # media_descr_list = j["medium_list"]
    return j


def writejson(data=None, path="data/config.json"):
    j = json.dumps(data)
    with open(path) as f:
        f.write(j)


# def read_transducer(json_str=None, json_path=None):

