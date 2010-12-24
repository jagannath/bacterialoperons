#! usr/bin/python










if __name__ == '__main__':
    
    [argument] = sys.argv[1:]
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - gen_fgoc_calc.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    