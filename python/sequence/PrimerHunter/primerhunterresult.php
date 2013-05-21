<?php
$max_filesize = ini_get('upload_max_filesize');
$targets_tmp = $_FILES['tf']['tmp_name'];
$nontargets_tmp = $_FILES['nf']['tmp_name'];
$primers_tmp = $_FILES['pf']['tmp_name'];
$mailto = $_POST["mailto"];
$id = $_GET["id"];
$processDone=0;
$errorMessage="";
if ($_FILES['tf']['error'] == UPLOAD_ERR_INI_SIZE)
{
	$errorMessage = "The size of the targets file is larger than the maximum allowed by the web server ($max_filesize)<BR>";
	$processDone = -2;
}
if ($_FILES['nf']['error'] == UPLOAD_ERR_INI_SIZE)
{
	$errorMessage = "The size of the non targets file is larger than the maximum allowed by the web server ($max_filesize)<BR>";
	$processDone = -2;
}
if ($_FILES['pf']['error'] == UPLOAD_ERR_INI_SIZE)
{
	$errorMessage = "The size of the primers file is larger than the maximum allowed by the web server ($max_filesize)<BR>";
	$processDone = -2;
}
if($processDone!=-2 && empty($id)&& empty($targets_tmp) && empty($primers_tmp)) 
{
	$errorMessage = "The server has detected an internal error. Try again later or contact the administrator.<BR>";
	$processDone = -2;
}
if($processDone!=-2 && (!empty($targets_tmp) || !empty($primers_tmp)))
{
	$ip = $_SERVER['REMOTE_ADDR'];
	$id="";
	$pos = strrpos($ip, ".");
	while($pos!==false)
	{
		$id = $id . substr($ip,$pos+1);
		$ip = substr($ip,0,$pos);
		$pos = strrpos($ip, ".");
	}
	$id = $id . $ip;
	$id = $id . date("YmdHisu");
	$cmd = "( nice ./primerhunter";
	if (!empty($targets_tmp)) {
		$targetsFile = "input/targets$id.txt";
		$cmd = $cmd ." -tf $targetsFile";
	}
	if (!empty($nontargets_tmp)) {
		$nonTargetsFile = "input/nonTargets$id.txt";
		$cmd = $cmd ." -nf $nonTargetsFile";
	}
	if (!empty($primers_tmp)) {
		$primersFile = "input/primers$id.txt";
		$cmd = $cmd ." -pf $primersFile";
	}
	$cmd = $cmd ." -pof output/primers$id.txt";

       	
	foreach ($_POST as $key=>$val) {
		if (!empty($val) && $key != "mailto") {
			$cmd = $cmd . " -$key ";
			if ($key != "full_stats") {
				$cmd = $cmd . $val;
			}
		}
	}
	$cmd = $cmd . " > output/output$id.txt";

	//Write mail info file
	if(!empty($mailto) && $mailto != "NONE")
	{
		$cmd = $cmd . "; ";
		$pos = strrpos($_SERVER ["REQUEST_URI"], "/");
		$outputURL = "http://".$_SERVER ["SERVER_NAME"].substr($_SERVER ["REQUEST_URI"],0,$pos)."/primerhunterresult.php?id=$id";
		$mailfile = fopen("output/mail$id.txt","wt");
		fwrite ($mailfile,"to: $mailto\n");
		fwrite ($mailfile,"subject: Primer Hunter Results\n");
		fwrite ($mailfile,"\nYour results are available in $outputURL\n");
		fclose ($mailfile);
		$cmd = $cmd . "/usr/sbin/sendmail -t < output/mail$id.txt  ";
	}
	$cmd = $cmd . ";touch output/done$id.txt) >& output/log$id.txt  &";
	//echo $cmd;	
	if (!empty($targets_tmp)) {
		move_uploaded_file($targets_tmp, $targetsFile);
	}
	if (!empty($nontargets_tmp)) {
		move_uploaded_file($nontargets_tmp, $nonTargetsFile);
	}
	if (!empty($primers_tmp)) {
		move_uploaded_file($primers_tmp, $primersFile);
	}
	system($cmd);
	$processDone=-1;
}
if($processDone>=0 && file_exists("output/done".$id.".txt"))
{
	$processDone = 1;
}

if($processDone==1) {
	$styleProcessing="display:none";
	$styleDone="display:";
	$styleError="display:none";
} else if ($processDone != -2) {
	$styleProcessing="display:";
	$styleDone="display:none";
	$styleError="display:none";
} else {
	$styleProcessing="display:none";
	$styleDone="display:none";
	$styleError="display:";

}

?>
<html>
  <head>
    <title>Primer Hunter</title>
<?
if($processDone==0)
{	
?>
    <meta http-equiv="Refresh" content="10; url=primerhunterresult.php?id=<?=$id?>">
<?
}
else if ($processDone == -1)
{
?>
    <meta http-equiv="Refresh" content="0; url=primerhunterresult.php?id=<?=$id?>"> 
<?
}
?>
  </head>
  <body>
    <table>
      <tr>
	<td align="center">
	  <h2>Primer Hunter</h2>
	</td>
      </tr>
      <tr style="<?=$styleDone?>">
        <td>
		Your request has been processed. See below your results. Primers output file is available <a href="output/primers<?=$id?>.txt">here</a>. You can load this file <a href=primerhunter.php?refilter=1>here</a> to refilter selected primers by changing some design parameters.               
        </td>
      </tr>
      <tr style="<?=$styleProcessing?>">
        <td>
		Your request have been received and it is being processed. If you provided an e-mail address you will receive an email as soon as the results of your request become available. Otherwise, you can bookmark this page to see your results later.
        </td>
      </tr>
      <tr style="<?=$styleError?>">
        <td>
		<?=$errorMessage?>
        </td>
      </tr>

      <tr>
        <td>
		<a href="primerhunter.php">Go back to design primers</a>
        </td>
      </tr>
      <tr style="<?=$styleDone?>">
        <td>
		<pre>
<?
if($processDone==1)
{
	$outputfile = fopen("output/output".$id.".txt","rt");
	while (!feof($outputfile)) {
		$line = fgets($outputfile);
		echo $line;
	}
	fclose($outputfile);
}
?>
		</pre>	
        </td>
      </tr>
  </body>
</html>

