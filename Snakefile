configfile: "config.yaml"
#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare wildcards and variables
data_dir = config["all"]["data_dir"]
output_dir = config["all"]["output_dir"]
#-------------------------------------------------------------------------------------------------------------------
# 0.2 specify target rules
rule all:
    input:
        output_dir + 'path'
        
#+++++++++++++++++++++++++++++++++++++++++ 1 XXXXXXXXXXXXXXXXXXXXXXXX  +++++++++++++++++++++++++++++++++++++++++++++
# 1.1  
rule RuleName:
    input:
        data_dir + 'path',
    output:
       output_dir + 'path'
    conda:
       "envs/R.yaml"
    script:
        "scripts/"

#-------------------------------------------------------------------------------------------------------------------
# 1.2
