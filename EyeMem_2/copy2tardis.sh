#!/bin/bash


## All IDs 

#sub-11009 sub-11012 sub-11018 sub-11019 sub-11021 sub-11022 sub-11024 sub-11029 sub-11031 sub-11032 sub-11038 sub-11046 sub-11048 sub-11049 sub-11051 sub-11060 
#sub-11067 sub-11074 sub-11075 sub-11083 sub-11092 sub-11102 sub-12005 sub-12014 sub-12016 sub-12030 sub-12033 sub-12034 sub-12036 sub-12037 sub-12039 sub-12050 
#sub-12053 sub-12054 sub-12070 sub-12084 sub-12090 sub-12091 sub-12100 sub-12106 sub-21006 sub-21008 sub-21010 sub-21026 sub-21035 sub-21044 sub-21047 sub-21052 
#sub-21059 sub-21065 sub-21069 sub-21071 sub-21072 sub-21079 sub-21081 sub-21086 sub-21088 sub-21093 sub-21094 sub-21095 sub-21104 sub-21105 sub-21114 sub-21115 
#sub-21116 sub-22004 sub-22007 sub-22011 sub-22020 sub-22021 sub-22023 sub-22025 sub-22040 sub-22041 sub-22043 sub-22045 sub-22056 sub-22058 sub-22061 sub-22062 
#sub-22064 sub-22066 sub-22068 sub-22076 sub-22080 sub-22099 sub-22103

server_dir="$HOME/eyemem/BIDS/B_data"
tardis_dir="$HOME/eyemem_tardis/BIDS/B_data"
# still to copy:
id_list="sub-21059 sub-21065 sub-21069 sub-21071 sub-21072 sub-21079 sub-21081 sub-21086 sub-21088 sub-21093 sub-21094 sub-21095 sub-21104 sub-21105 sub-21114 sub-21115 sub-21116 sub-22004 sub-22007 sub-22011 sub-22020 sub-22021 sub-22023 sub-22025 sub-22040 sub-22041 sub-22043 sub-22045 sub-22056 sub-22058 sub-22061 sub-22062 sub-22064 sub-22066 sub-22068 sub-22076 sub-22080 sub-22099 sub-22103";

for id in $id_list; do 
    data=$server_dir/$id
    echo "Copying $id"
    cp -r $data $tardis_dir 
done