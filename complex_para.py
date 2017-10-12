def main():
    import argparse
    from tools_para import message,load_model,analyze_results
    #import time
    #from time import *
    import timeit
    print
    print
    nprocess = 1
    nn = 1e3
    parser =  argparse.ArgumentParser(description='Options to run the antidune model')
    parser.add_argument('-m','--mode', help='Mode: d-Development,p-Production',required=True)
    parser.add_argument('-s','--sample',help='M-C Sample size',required=False)
    parser.add_argument('-p','--processors',help='# number of processors',required=False)
    args = parser.parse_args()

    if str(args.processors) != 'None':
        nprocess=int(args.processors)

    if str(args.sample) != 'None':
        nn=int(args.sample)

    mc_n = int(nn)
    if str(args.mode) == 'd':
        from progress.bar import Bar
        message()
        mode=str(args.mode)
        print "Number of processors: \t\t\t{t1:3d}".format(t1=nprocess)
        print "Sample size for M-C simulations:  {t2:9d}".format(t2=mc_n)
    if str(args.mode) == 'p':
        mode = str(args.mode)
    start_time = timeit.default_timer()
    data=load_model(mc_n,nprocess)
    analyze_results(data)
    elapsed = timeit.default_timer() - start_time
    print('Run time = {test1:6.1f} seconds\n'.format(test1=elapsed))
    # TODO: Implement choice as to whether density should be calculated or read rom file
    # FIXME: afjdkajfa fjljk
    #plot_data(save_fig=0,plot_density=1)
    print
    print "\t\t     \033[1m D O N E !\033[0m"



if __name__ == '__main__':
	print "changed some thing"
	print "even more"
    main()
