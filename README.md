# calmet_dump
dump data from calmet format 3d met data

## prerequisite

Need fortran compiler.  Intel and gfortran are known to work on Linux

## compile
`make`
This creates executable calmet_dump

## usage

### 1. simple dump

`calmet_dump calmet.m3d`
`calmet_dump 2, 3 calmet.m3d`

The above dumps all variables, all tstep for grid cell (i,j) = (1,1), or (i,j) = (2,3), for example

### 2. a bit more formatted output

`calmet_dump -`

The above commane ('-' as argument) will prompt user for further inputs.  it will ask
1. filename
2. idx jdx to extract (default = 1,1)
3. start time to extract, YYYYJJJHH, SSSS (default = 1999010100, 0000, i.e. beg of year 1900)
4. stop time to extract, YYYYJJJHH, SSSS (default = 2099123123, 3599, i.e. end of year 2099)
5. list of values to extract, one value in each line, end with blank line

To got time and variable in the file, i first run with simple run method, write donw what you want to do, and then type (or you can write it down into a text file and redirect into code like `calmet_dump - < my_inputs > out.csv`) 

Example goes  like below.   I redirected output to go to out.csv.  code asks me each input and i have to type each without any error.  then i got csv output.

```
$ ./calmet_dump - > out.csv
fname: calmetv6.freestone.20190910_20190919.m3d
 calmetv6.freestone.20190910_20190919.m3d
i,j: 2, 3
           2           3
btime (yyyyjjjhh ssss): 201926005 0
   201926005           0
etime (yyyyjjjhh ssss): 201926020 0
   201926020           0
variables (blank to all): U-LEV
 more variable (blank to end): V-LEV
 more variable (blank to end): USTAR
 more variable (blank to end): IPGT
 more variable (blank to end):
 finished, n=           4
 U-LEV
 V-LEV
 USTAR
 IPGT
$ head out.csv
U-LEV   ,201926105,   0,201926106,   0,   2,   3,   1,-0.7965069E+00
V-LEV   ,201926105,   0,201926106,   0,   2,   3,   1,-0.1288582E+01
IPGT    ,201926105,   0,201926106,   0,   2,   3,   1,             4
USTAR   ,201926105,   0,201926106,   0,   2,   3,   1, 0.6000000E-01
U-LEV   ,201926106,   0,201926107,   0,   2,   3,   1,-0.2347563E+01
V-LEV   ,201926106,   0,201926107,   0,   2,   3,   1,-0.2007854E+01
IPGT    ,201926106,   0,201926107,   0,   2,   3,   1,             4
USTAR   ,201926106,   0,201926107,   0,   2,   3,   1, 0.1600000E+00
U-LEV   ,201926107,   0,201926108,   0,   2,   3,   1,-0.2188389E+01
V-LEV   ,201926107,   0,201926108,   0,   2,   3,   1,-0.2298859E+01
...
$ 
```

### 3. netCDF format

The command below creates `m3dfile.nc` which has most of data from m3d fie

`calmet_dump -n m3dfile.m3d`

or below, if you want to use different name (or location)

`calmet_dump -n m3dfile.m3d name_that_i_want.nc`

## Next

netcdf file lacks dataset level attributes, such as projection information, grid information etc.  I may get back to it if i need them for some reason...
