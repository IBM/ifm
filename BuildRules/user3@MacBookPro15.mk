# For macOS Catalina (v. 10.15.5)
USER_LDFLAGS  = -lomp
OPTFLAG       = -g -O3 -Xpreprocessor -fopenmp # -Xpreprocessor must precede to -fopenmp
USER_CXXFLAGS = -Wno-format # To avoid tons of warning messages
