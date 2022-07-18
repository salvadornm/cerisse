case  $1 in
        backup )
        tput setaf 2
        echo " backing up.. in $2 "
	mkdir $2
	cp screen $2
	cp fort.* $2
	cp inputs $2	
	cp *.F90 $2
	echo " .....   [done] "
	tput sgr0
    ;;
	restart)
	echo " prepare restart from step $2  [NRY]"
	tput setaf 2	
	echo " .....   [done] "
	tput sgr0
	;;
	plot)
	echo " creating AMR plots ..."
	ls -1 plt*/Header | tee movie.visit
	;;
	kill )
	killall CNS*
	;;
	check )
	ps aux | grep compreal | grep -v grep
	;;
	clean )
	rm -rf plt*
	rm -f screen
	rm -f *.visit
	rm -f Backtra*
	;;
	purge )
	rm -rf plt*
        rm -f screen
        rm -f *.visit
        rm -f Backtra*
	make clean
	;;
    *)
	tput setaf 2
	echo " CNS running .."
	tput sgr0
        ps aux | grep CNS* | grep -v grep
esac
