import netgen.config

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--cmake-dir", help="print path to CMake config files", action='store_true')
    args = parser.parse_args()
    if(args.cmake_dir):
        print(netgen.config.get_cmake_dir())

