import gzip

import argparse


try:
    parser = argparse.ArgumentParser(description='Please give the name of the input file.')
    parser.add_argument("--ncores", help="Number of cores", type=int,
                        action="store")
    parser.add_argument('inputfile',
                        metavar='File', type=str, help='Input file name')
    args = parser.parse_args()
except:
    print("Please give an input file")
    raise SystemExit

if args.ncores:
    n_cores=args.ncores
else:
    import multiprocessing as mp
    n_cores=mp.cpu_count()

#bannerlines=[]
#eventlines=[]


infile=args.inputfile




## infile should be something.lhe.gz

fhs=[]

try:
    for x in range(n_cores):
        os.remove('Split/split_'+str(x)+'.lhe.gz')
except:
    pass

try:
    for x in range(n_cores):
        tfh=gzip.open('Split/split_'+str(x)+'.lhe.gz','wb')
        fhs.append(tfh)
except:
    print("Could not open output files")
    raise SystemExit

try:
    with gzip.open(infile, 'rb') as f:
        donebanner=False
        done=False
        eventhandle=0
        for line in f:
            sline=line.decode('ascii')
            if not donebanner:
                if '<event' in sline:
                    donebanner=True
                    print('Found first event line')
                    print(sline.strip())
                    fhs[eventhandle].write(line)
#                eventlines=[line]

                else:
                    print(sline.strip())
                    for fh in fhs: ## now write the banner to all of the output files
                        fh.write(line)
            else:
                if not done:
                    if '</event>' in sline:
                        # finish old event and advance the counter
                        fhs[eventhandle].write(line)
                        eventhandle=(eventhandle+1) % n_cores
                    else:
                    
                        if '</LesHouchesEvent' in sline:
                            done=True
                            for fh in fhs:
                                fh.write(line)
                                fh.close()
                        else:
                            # regular line
                            fhs[eventhandle].write(line)
                
except:
    print("Could not open input files")
    raise SystemExit
