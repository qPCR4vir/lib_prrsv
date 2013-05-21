#!/bin/csh
foreach x ( $argv )
    tail -1 $x | grep -q '^Mega BLAST run finished, ' 
    if( $status != 0 ) then
	set t = `echo $x | sed 's/.*\.//'`
	echo $t
	#rm $x
    endif
end
