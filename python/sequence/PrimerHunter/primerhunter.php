<?php
$max_filesize = ini_get('upload_max_filesize');
//Read constants from primerhunter.h
$constantsfile = fopen("primerhunter.h","rt");
while (!feof($constantsfile)) {
	$line = fgets($constantsfile);
	$tok = strtok($line, " \n\t");
	if ($tok !== false && $tok == "#define") {
		$ckey = strtok(" \n\t");
		$cvalue = strtok(" \n\t");
		$primerconstants [$ckey] = $cvalue;
	}   
}
fclose($constantsfile);
if($_GET["refilter"]==1) {
	$refilterJavascript = "true";
	$styleFirst="display:none";
	$styleRefilter="display:";
} else {
	$refilterJavascript = "false";
	$styleFirst="display:";
	$styleRefilter="display:none";
}
?>

<html>
  <head>
    <title>Primer Hunter</title>
    <script language="javascript">
	function validateForm (form) {
		refilter = <?=$refilterJavascript?>;
		message = "";
		if (!refilter && !form.tf.value )  {
			message += "Targets file required\n";
		}
		if(refilter && !form.pf.value) {
			message += "Primers file required\n";
		}
		if (!refilter && !form.minPrimerLength.value) {
			message += "Min Primer Length required\n";
		} else if (!refilter && !form.minPrimerLength.value.match(/^([0-9]+)$/)) {
			message += "Min Primer Length must be numeric\n";
		}
		if (!refilter && !form.maxPrimerLength.value) {
			message += "Max Primer Length required\n";
		} else if (!refilter && !form.maxPrimerLength.value.match(/^([0-9]+)$/)) {
			message += "Max Primer Length must be numeric\n";
		}
		if (!form.minProdLength.value) {
			message += "Min Product Length required\n";
		} else if (!form.minProdLength.value.match(/^([0-9]+)$/)) {
			message += "Min Product Length must be numeric\n";
		}
		if (!form.maxProdLength.value) {
			message += "Max Product Length required\n";
		} else if (!form.maxProdLength.value.match(/^([0-9]+)$/)) {
			message += "Max Product Length must be numeric\n";
		}
		if (!refilter && !form.forwardPrimer.value) {
			message += "Forward Primer required\n";
		} else if (!refilter && form.forwardPrimer.value!="NONE" && !form.forwardPrimer.value.match(/^([GATCN]+)$/)) {
			message += "Forward Primer must be NONE or have just G, A, T, C or N as characters\n";
		}
		if (!refilter && !form.reversePrimer.value) {
			message += "Reverse Primer required\n";
		} else if (!refilter && form.reversePrimer.value!="NONE" && !form.reversePrimer.value.match(/^([GATCN]+)$/)) {
			message += "Reverse Primer must be NONE or have just G, A, T, C or N as characters\n";
		}
		if (!refilter && !form.beginPosForward.value) {
			message += "Begin of range of positions for forward primers required\n";
		} else if (!refilter && !form.beginPosForward.value.match(/^([0-9]+)$/)) {
			message += "Begin of range of positions for forward primers must be numeric\n";
		}
		if (!refilter && !form.endPosForward.value) {
			message += "End of range of positions for forward primers required\n";
		} else if (!refilter && !form.endPosForward.value.match(/^([0-9]+)$/)) {
			message += "End of range of positions for forward primers must be numeric\n";
		}
		if (!refilter && !form.beginPosReverse.value) {
			message += "Begin of range of positions for reverse primers required\n";
		} else if (!refilter && !form.beginPosReverse.value.match(/^([0-9]+)$/)) {
			message += "Begin of range of positions for reverse primers must be numeric\n";
		}
		if (!refilter && !form.endPosReverse.value) {
			message += "End of range of positions for reverse primers required\n";
		} else if (!refilter && !form.endPosReverse.value.match(/^([0-9]+)$/)) {
			message += "End of range of positions for reverse primers must be numeric\n";
		}
		if (!refilter && !form.tmask.value) {
			message += "Targets Mask required\n";
		} else if (!refilter && !form.tmask.value.match(/^([01]+)$/)) {
			message += "Targets Mask must be binary\n";
		}
		if (!refilter && !form.nmask.value) {
			message += "Non Targets Mask required\n";
		} else if (!refilter && form.nmask.value!="NONE" && !form.nmask.value.match(/^([01]+)$/)) {
			message += "Non Targets Mask must be binary or NONE\n";
		}
		if (!refilter && !form.dmask.value) {
			message += "Degeneracy Mask required\n";
		} else if (!refilter && !form.dmask.value.match(/^([14]+)$/)) {
			message += "Degeneracy Mask must have just 1 or 4 as digits\n";
		}
		if (!refilter && !form.numSourceSeq.value) {
			message += "Number of Source Sequences required\n";
		} else if (!refilter && !form.numSourceSeq.value.match(/^([0-9]+)$/)) {
			message += "Number of Source Sequences must be numeric\n";
		}
		if (!form.minCoverageTargets.value) {
			message += "Minimum Targets Coverage required\n";
		} else if (!form.minCoverageTargets.value.match(/(^-?\d\d*\.\d*$)|(^-?\d\d*$)|(^-?\.\d\d*$)/)) {
			message += "Minimum Targets Coverage must be numeric\n";
		} 
		if (!form.maxCoverageNonTargets.value) {
			message += "Maximum Non Targets Coverage required\n";
		} else if (!form.maxCoverageNonTargets.value.match(/(^-?\d\d*\.\d*$)|(^-?\d\d*$)|(^-?\.\d\d*$)/)) {
			message += "Maximum Non Targets Coverage must be numeric\n";
		}
		if (!form.maxSelfScore.value) {
			message += "Max self score required\n";
		} else if (!form.maxSelfScore.value.match(/^([0-9]+)$/)) {
			message += "Max self score must be numeric\n";
		} 
		if (!form.maxEndScore.value) {
			message += "Max end score required\n";
		} else if (!form.maxEndScore.value.match(/^([0-9]+)$/)) {
			message += "Max end score must be numeric\n";
		}
		if (!form.minGCContent.value) {
			message += "Min GC Content required\n";
		} else if (!form.minGCContent.value.match(/(^-?\d\d*\.\d*$)|(^-?\d\d*$)|(^-?\.\d\d*$)/)) {
			message += "Min GC Content must be numeric\n";
		} 
		if (!form.maxGCContent.value) {
			message += "Max GC Content required\n";
		} else if (!form.maxGCContent.value.match(/(^-?\d\d*\.\d*$)|(^-?\d\d*$)|(^-?\.\d\d*$)/)) {
			message += "Max GC Content must be numeric\n";
		}
		if (!form.gcClamp.value) {
			message += "GC Clamp required\n";
		} else if (!form.gcClamp.value.match(/^([0-9]+)$/)) {
			message += "GC Clamp must be numeric\n";
		}
		if (!form.maxPolyX.value) {
			message += "Max Poly-X required\n";
		} else if (!form.maxPolyX.value.match(/^([0-9]+)$/)) {
			message += "Max Poly-X must be numeric\n";
		}
		if (!form.primersConc.value) {
			message += "Primer Concentration required\n";
		} else if (!form.primersConc.value.match(/(^-?\d\d*\.\d*$)|(^-?\d\d*$)|(^-?\.\d\d*$)/)) {
			message += "Primer Concentration must be numeric\n";
		} 
		if (!form.templateConc.value) {
			message += "Template Concentration required\n";
		} else if (!form.templateConc.value.match(/(^-?\d\d*\.\d*$)|(^-?\d\d*$)|(^-?\.\d\d*$)/)) {
			message += "Template Concentration must be numeric\n";
		}
		if (!form.saltConc.value) {
			message += "Salt Concentration required\n";
		} else if (!form.saltConc.value.match(/(^-?\d\d*\.\d*$)|(^-?\d\d*$)|(^-?\.\d\d*$)/)) {
			message += "Salt Concentration must be numeric\n";
		}   
		if (!form.minTempTargets.value) {
			message += "Targets Min Melting Temp required\n";
		} else if (!form.minTempTargets.value.match(/(^-?\d\d*\.\d*$)|(^-?\d\d*$)|(^-?\.\d\d*$)/)) {
			message += "Targets Min Melting Temp must be numeric\n";
		} 
		if (!form.maxTempTargets.value) {
			message += "Targets Max Melting Temp required\n";
		} else if (!form.maxTempTargets.value.match(/(^-?\d\d*\.\d*$)|(^-?\d\d*$)|(^-?\.\d\d*$)/)) {
			message += "Targets Max Melting Temp must be numeric\n";
		}
		if (!form.maxTempNonTargets.value) {
			message += "Non Targets Max Melting Temp required\n";
		} else if (!form.maxTempNonTargets.value.match(/(^-?\d\d*\.\d*$)|(^-?\d\d*$)|(^-?\.\d\d*$)/)) {
			message += "Non Targets Max Melting Temp must be numeric\n";
		} 
		if (!form.deltaTempNonTargets.value) {
			message += "Non Targets Delta Melting Temp required\n";
		} else if (!form.deltaTempNonTargets.value.match(/(^-?\d\d*\.\d*$)|(^-?\d\d*$)|(^-?\.\d\d*$)/)) {
			message += "Non Targets Delta Melting Temp must be numeric\n";
		}
		if (!form.maxPairTempDiff.value) {
			message += "Max Pair Tm Diff required\n";
		} else if (!form.maxPairTempDiff.value.match(/(^-?\d\d*\.\d*$)|(^-?\d\d*$)|(^-?\.\d\d*$)/)) {
			message += "Max Pair Tm Diff must be numeric\n";
		} 
		if (!form.primersLabel.value) {
			message += "Primers Label required\n";
		}
		if(message.length >0) {
			alert(message);
			return false;
		}
		return true;
	}	
    </script>	
  </head>
  <body>
  <form enctype="multipart/form-data" action="primerhunterresult.php" method="post" onsubmit="return validateForm(this);">
    <table>
      <tr>
	<td align="center" colspan="4">
	  <h2>Primer Hunter</h2>
	</td>
      </tr>
      <tr style="<?=$styleFirst?>">
        <td colspan="4">
                <b>New Users:</b> Load your targets file and (optional) non targets file in fasta format in the corresponding boxes and click on the design primers button.
        </td>
      </tr>
      <tr style="<?=$styleRefilter?>">
        <td colspan="4">
                <b>New Users:</b> Load your primers summary file in the primers file box and click on the design primers button.
        </td>
      </tr>
      <tr>
        <td colspan="4">
		Optionally, if you want to be notified as soon as the process is done, write at the end an email address where you want the results to be sent.
        </td>
      </tr>

      <tr>
        <td colspan="4">
                Read the <a href="phdocs/index.html">documentation</a> for further information.
        </td>
      </tr>
      <tr style="<?=$styleFirst?>">
        <td colspan="4">
                <a href="primerhunter.php?refilter=1">Go to refilter primers</a> <a href="phdocs/index.html#twostepsdesign" style="font-size:70%">What is this?</a>
        </td>
      </tr>
      <tr style="<?=$styleRefilter?>">
        <td colspan="4">
                <a href="primerhunter.php">Go back to design primers</a>
        </td>
      </tr>

      <tr style="<?=$styleFirst?>">
	<td><a href="phdocs/index.html#tf">Targets file</a>(Max <?=$max_filesize?>)
	</td>
	<td colspan="2"><input type="file" name="tf" /> 
	</td>
	<td><a href="sampleTarget2.txt">Sample</a>
	</td>
      </tr>
      <tr style="<?=$styleFirst?>">
	<td><a href="phdocs/index.html#nf">Non Targets file</a>(Max <?=$max_filesize?>)
	</td>
	<td colspan="2"><input type="file" name="nf" /> 
	</td>
	<td><a href="sampleNonTarget2.txt">Sample</a>
	</td>
      </tr>
      <tr style="<?=$styleRefilter?>">
	<td><a href="phdocs/index.html#pf">Primers file</a>(Max <?=$max_filesize?>)
	</td>
	<td colspan="2"><input type="file" name="pf" /> 
	</td>
	<td><a href="samplePrimer2.txt">Sample</a>
	</td>
      </tr>

      <tr style="<?=$styleFirst?>">
	<td><a href="phdocs/index.html#minPrimerLength">Min Primer Length</a>
	</td>
	<td><input type="text" name="minPrimerLength" value="<?echo $primerconstants ['DEF_MIN_PRIMER_LENGTH']?>" /> 
	</td>
	<td><a href="phdocs/index.html#maxPrimerLength">Max Primer Length</a>
	</td>
	<td><input type="text" name="maxPrimerLength" value="<?echo $primerconstants ['DEF_MAX_PRIMER_LENGTH']?>"/>
	</td>
      </tr>
      
      <tr>
	<td><a href="phdocs/index.html#minProdLength">Min Product Length</a>
	</td>
	<td><input type="text" name="minProdLength" value="<?echo $primerconstants ['DEF_MIN_PROD_LENGTH']?>"/> 
	</td>
	<td><a href="phdocs/index.html#maxProdLength">Max Product Length</a>
	</td>
	<td><input type="text" name="maxProdLength" value="<?echo $primerconstants ['DEF_MAX_PROD_LENGTH']?>"/>
	</td>
      </tr>
      <tr style="<?=$styleFirst?>">
	<td><a href="phdocs/index.html#forwardPrimer">Forward Primer (5'-3')</a>
	</td>
	<td><input type="text" name="forwardPrimer" value=<?echo $primerconstants ['DEF_FORWARD_PRIMER']?>/> 
	</td>
	<td><a href="phdocs/index.html#reversePrimer">Reverse Primer (5'-3')</a>
	</td>
	<td><input type="text" name="reversePrimer" value=<?echo $primerconstants ['DEF_REVERSE_PRIMER']?>/> 
	</td>
      </tr>
      <tr style="<?=$styleFirst?>">
	<td colspan=4><a href="phdocs/index.html#beginPosForward">Range of positions for forward primers</a>
	</td>
      </tr>
      <tr style="<?=$styleFirst?>">
	<td>Begin:
	</td>
	<td><input type="text" name="beginPosForward" value="<?echo $primerconstants ['DEF_BEGINPOS_FORWARD']?>"/> 
	</td>
	<td>End:
	</td>
	<td><input type="text" name="endPosForward" value="<?echo $primerconstants ['SEQUENCES_MAX_SIZE']?>"/>
	</td>
      </tr>
      <tr style="<?=$styleFirst?>">
	<td colspan=4><a href="phdocs/index.html#beginPosReverse">Range of positions for reverse primers</a>
	</td>
      </tr>
      <tr style="<?=$styleFirst?>">
	<td>Begin:
	</td>
	<td><input type="text" name="beginPosReverse" value="<?echo $primerconstants ['DEF_BEGINPOS_REVERSE']?>"/> 
	</td>
	<td>End:
	</td>
	<td><input type="text" name="endPosReverse" value="<?echo $primerconstants ['SEQUENCES_MAX_SIZE']?>"/>
	</td>
      </tr>
      <tr style="<?=$styleFirst?>">
	<td><a href="phdocs/index.html#tmask">Targets Mask (3'-5')</a>
	</td>
	<td><input type="text" name="tmask" value=<?echo $primerconstants ['DEF_MASK_TARGETS']?>/> 
	</td>
	<td><a href="phdocs/index.html#nmask">Non Targets Mask (3'-5')</a>
	</td>
	<td><input type="text" name="nmask" value=<?echo $primerconstants ['DEF_MASK_NONTARGETS']?>/> 
	</td>
      </tr>
      <tr style="<?=$styleFirst?>">
	<td><a href="phdocs/index.html#dmask">Degeneracy Mask (3'-5')</a>
	</td>
	<td><input type="text" name="dmask" value=<?echo $primerconstants ['DEF_DEGENERACY']?>/>
	</td>
	<td><a href="phdocs/index.html#numSourceSeq">Number of Source Sequences</a>
	</td>
	<td><input type="text" name="numSourceSeq" value="<?echo $primerconstants ['DEF_SEQ_FOR_CANDS']?>"/>
	</td>
      </tr>
      <tr>
	<td><a href="phdocs/index.html#minCoverageTargets">Minimum Targets Coverage (%)</a>
	</td>
	<td><input type="text" name="minCoverageTargets" value="<?echo $primerconstants ['DEF_MIN_COVERAGE_TARGETS']?>"/>
	</td>
	<td><a href="phdocs/index.html#maxCoverageNonTargets">Maximum Non Targets Coverage (%)</a>
	</td>
	<td><input type="text" name="maxCoverageNonTargets" value="<?echo $primerconstants ['DEF_MAX_COVERAGE_NONTARGETS']?>"/>
	</td>
      </tr>

      <tr>
	<td><a href="phdocs/index.html#maxSelfScore">Max self score</a>
	</td>
	<td><input type="text" name="maxSelfScore" value="<?echo $primerconstants ['DEF_MAX_SELF_SCORE']?>"/> 
	</td>
	<td><a href="phdocs/index.html#maxEndScore">Max end score</a>
	</td>
	<td><input type="text" name="maxEndScore" value="<?echo $primerconstants ['DEF_MAX_END_SCORE']?>"/> 
	</td>
      </tr>
      <tr>
	<td><a href="phdocs/index.html#minGCContent">Min GC Content</a>
	</td>
	<td><input type="text" name="minGCContent" value="<?echo $primerconstants ['DEF_MIN_GCCONTENT']?>"/> 
	</td>
	<td><a href="phdocs/index.html#maxGCContent">Max GC Content</a>
	</td>
	<td><input type="text" name="maxGCContent" value="<?echo $primerconstants ['DEF_MAX_GCCONTENT']?>"/> 
	</td>
      </tr>
      <tr>
	<td><a href="phdocs/index.html#gcClamp">GC Clamp</a>
	</td>
	<td><input type="text" name="gcClamp" value="<?echo $primerconstants ['DEF_GCCLAMP']?>"/> 
	</td>
	<td><a href="phdocs/index.html#maxPolyX">Max Poly-X</a>	
	</td>
	<td><input type="text" name="maxPolyX" value="<?echo $primerconstants ['DEF_MAX_POLY_X']?>"/> 
	</td>
      </tr>

      <tr>
	<td><a href="phdocs/index.html#primersConc">Primer concentration (M)</a>
	</td>
	<td><input type="text" name="primersConc" value="<?echo $primerconstants ['DEF_CONC_PRIMERS']?>"/> 
	</td>
	<td><a href="phdocs/index.html#templateConc">Template Concentration (M)</a>
	</td>
	<td><input type="text" name="templateConc" value="<?echo $primerconstants ['DEF_CONC_SEQUENCES']?>"/> 
	</td>
      </tr>
      <tr>
	<td><a href="phdocs/index.html#saltConc">Salt Concentration (M)</a>
	</td>
	<td><input type="text" name="saltConc" value="<?echo $primerconstants ['DEF_SALT']?>"/>
	</td>
	<td><a href="phdocs/index.html#saltCorrMethod">Salt Correction Method</a>
	</td>
	<td><select name="saltCorrMethod" ><option value=1 selected>Santalucia</option><option value=2 >Owczarzy</option></select>
	</td>
      </tr>
      <tr>
	<td><a href="phdocs/index.html#minTempTargets">Targets Min Melting Temp (C)</a>
	</td>
	<td><input type="text" name="minTempTargets" value="<?echo $primerconstants ['DEF_MIN_TEMP_TARGETS']?>"/>
	</td>
	<td><a href="phdocs/index.html#maxTempTargets">Targets Max Melting Temp (C)</a>
	</td>
	<td><input type="text" name="maxTempTargets" value="<?echo $primerconstants ['DEF_MAX_TEMP_TARGETS']?>"/>
	</td>
      </tr>
      <tr>
	<td><a href="phdocs/index.html#maxTempNonTargets">Non Targets Max Melting Temp (C)</a>
	</td>
	<td><input type="text" name="maxTempNonTargets" value="<?echo $primerconstants ['DEF_MAX_TEMP_NONTARGETS']?>"/>
	</td>
	<td><a href="phdocs/index.html#deltaTempNonTargets">Non Targets Delta Melting Temp (C)</a>
	</td>
	<td><input type="text" name="deltaTempNonTargets" value="<?echo $primerconstants ['DEF_DELTA_TEMP_NONTARGETS']?>"/>
	</td>
      </tr>
      <tr>
	<td><a href="phdocs/index.html#maxPairTempDiff">Max Pair Tm Diff (C)</a>
	</td>
	<td><input type="text" name="maxPairTempDiff" value="<?echo $primerconstants ['DEF_MAX_PAIR_TEMP_DIFF']?>"/>	
	</td>

	<td>	
	</td>
	<td>
	</td>
      </tr>

      <tr>
	
	<td><a href="phdocs/index.html#primersLabel">Primers Label</a>
	</td>
	<td><input type="text" name="primersLabel" value=<?echo $primerconstants ['DEF_PRIMERS_LABEL']?>/>
	</td>
	<td><a href="phdocs/index.html#full_stats">Full stats</a>
	</td>
	<td><input type="checkbox" name="full_stats" />
	</td>

      </tr>

      <tr>
	<td><a href="phdocs/index.html#mailto">E-mail address</a>	
	</td>
	<td><input type="text" name="mailto" value="NONE"  />
	</td>
	<td>	
	</td>
	<td>
	</td>
      </tr>
      <tr>
	<td align="center" colspan="2">
	  <input type="submit" value="Design Primers"/>
	</td>
	<td align="center" colspan="2">
	  <input type="reset" value="Restore Defaults"/>
	</td>
      </tr>			
    </table>
  </form>
  </body>
</html>

