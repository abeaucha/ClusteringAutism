import subprocess

def execute_R(script, args):
    
    """
    Execute an R script.
    
    Arguments
    ---------
    script: str
        String containing the name of the R script to execute.
    args: dict
        Dictionary of key-value pairs containing command line arguments to pass to the script
    
    Returns
    -------
    None
    """
    
    args = [['--'+str(key), str(val)] for key, val in args.items()]
    args = sum(args, [])
    cmd = ['Rscript']+[script]+args
    subprocess.run(cmd)
    return