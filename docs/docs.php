<?php



include("./_header.php"); 

$path=getcwd()."/documents/";

$allpdfs=scandir($path);

unset($allpdfs[0]);

unset($allpdfs[1]);

echo "<ul>";

foreach ($allpdfs as $doc) :?>

	<li><a href="<?= $path.$doc;?>" target="_blank"><?= $doc;?></a></li>

<?php endforeach;

echo "</ul>";



include("./_footer.php"); 

?>