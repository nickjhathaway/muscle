ofiles=`echo *.o`

if [ ! $ofiles == "*.o" ] ; then
	rm *.o
fi

make -f Makefile > make.out 2> make.err

if [ ! $ofiles == "*.o" ] ; then
	rm *.o
fi

cat make.err
ls -l muscle
