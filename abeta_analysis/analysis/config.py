import os

# System settings
cluster = False

# munging settings
isomer_list = ["chiro"]
ratio_list = [15]

data_source = os.environ['PWD']
if cluster:
    data_source = '/dev/shm'