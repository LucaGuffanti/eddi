# A simple script to download vdw radii from the web

import requests
import pandas as pd
import regex as re
import subprocess
import os

elements = []

def download_periodic_table(verbose=False):
    if not os.path.exists('scripts/shielding/utils/periodic_table.csv'):
        if verbose:
            print("Downloading periodic table data")
        command = '''
        pwd
        cd utils/ &&
        wget https://gist.githubusercontent.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee/raw/1d92663004489a5b6926e944c1b3d9ec5c40900e/Periodic%2520Table%2520of%2520Elements.csv &&
        mv "Periodic Table of Elements".csv periodic_table.csv
        sed -i 's/Wolfram/Tungsten/g' periodic_table.csv
        sed -i 's/Cesium/Caesium/g' periodic_table.csv
        '''

        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        if process.returncode == 0 and verbose:
            print(f"Command executed successfully: {stdout.decode().strip()}")
        elif verbose:
            print(f"Error executing command: {stderr.decode().strip()}")


def get_data(verbose=False):
   
    download_periodic_table(verbose)
    
    df = pd.read_csv('utils/periodic_table.csv')[['AtomicNumber','Element', 'Symbol']]

    pattern = re.compile(r'van der Waals radius</a>: \d+')
    alt_pattern = re.compile(r'van der Waals radius</a>: \[ \d+ \]')

    length = len(df)

    for i in range(length):
        name = df['Element'][i].lower()
        if verbose:
            print(name)

        url = f"https://www.webelements.com/{name}/index.html"
        req = requests.get(url)

        if req.status_code == 200 and verbose:
            print(f"Page for {name} (Atomic Number: {df['AtomicNumber'][i]}) was found.")
        elif verbose:
            print(f"[ERROR] Page for {name} (Atomic Number: {df['AtomicNumber'][i]}) was not found")
            continue

        req_data = req.text
        match = pattern.search(req_data)
        if match and verbose:
            print(f"Match found for {name}: {match.group(0)}")
        elif verbose:
            print(f"Trying alternative pattern")
            match = alt_pattern.search(req_data)
            if match and verbose:
                print(f"Match found for {name}: {match.group(0)}")
            elif verbose:
                print(f"[ERROR] No match found for {name}")
                break

        radius = int(re.search(r'\d+', match.group(0)).group(0))
        elements.append({
            'AtomicNumber': df['AtomicNumber'][i], 
            'Element': df['Element'][i], 
            'Symbol': df['Symbol'][i], 
            'vdwRadius': radius
        })
        if verbose:
            print(f"===>> Radius: {radius}")

    output_df = pd.DataFrame(elements)
    try:
        os.mkdir("data/")
    except FileExistsError:
        pass

    output_df.to_csv('data/vdw_radii.csv', index=False)
            

def load_vdw_radii(verbose=False):
    if not os.path.exists('data/vdw_radii.csv'):
        get_data(verbose)

    df = pd.read_csv('data/vdw_radii.csv')
    return df['vdwRadius'].to_list()

def get_radii(verbose=False):

    download_periodic_table()

    # TODO: find the website containing all the nuclear radii
    # In order to parse the website response

    # Build a http request

    # Apply the necessary regex

    # Store data in a csv

def load_nuclear_radii(verbose=False):
    if not os.path.exists('data/nuclear_radii.csv'):
        get_radii(verbose)
    
    df = pd.read_csv('data/nuclear_radii.csv')
    return df['nuclearRadius']

if __name__ == "__main__":
    get_data(verbose=True)
