



def main():
    from tools import plot_data
    print
    print
    print "\t \t \033[1m Creating Plot for Antidunes \033[0m"
    print
    print
    # TODO: Implement choice as to whether density should be calculated or read rom file
    plot_data(save_fig=0,plot_density=1)
    print "\t\t\t \033[1m D O N E !\033[0m"



if __name__ == '__main__':
    main()
