#!/bin/csh
foreach x ( $argv )

    echo $x | grep -q '.dG$'
    if( $status != 0 ) then
	#set dgLsN = `ls -1 ${x}* | grep -c '.dG$'`
 	if (! -e ${x}.dG) then
	    echo $x
	endif
    endif
    

#    set dgCt = `grep -ch '^' $x`
#    set mbName = `echo $x | sed 's/.dG//' `
#    set mbCt = `grep -ch '^#' $mbName`
    
#    echo $dgCt $mbCt

#    if( $mbCt != $dgCt ) then
#	set t = `echo $mbName | sed 's/.*\.//'`
#	echo $t
#    else
#	echo $mbName
#    endif
#    echo
end
