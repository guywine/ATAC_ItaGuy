import pandas as pd
import gzip
import xml.etree.ElementTree as ET


if __name__=='__main__':
    with gzip.open('tables/Variation.xml.gz') as f_name:
        content = f_name.read()  
        tree = ET.parse(content)
        root = tree.getroot()
