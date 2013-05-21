#!/bin/csh
foreach x ( $argv )
    set dgCt = `grep -ch '^' $x`
    set mbName = `echo $x | sed 's/.dG//' `
    set mbCt = `grep -ch '^#' $mbName`
    
#    echo $dgCt $mbCt

    if( $mbCt != $dgCt ) then
	set t = `echo $mbName | sed 's/.*\.//'`
	echo $t
#    else
#	echo $mbName
    endif
#    echo
end
