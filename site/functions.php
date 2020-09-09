<?php



include("./_header.php"); 



?>



<div id='content'>



	<?php 

		echo "<ul>";

		$path=getcwd()."/functions";

		$funcs=scandir($path);



		unset($funcs[0]);

		unset($funcs[1]);

		foreach ($funcs as $func):

	?>

			<li>

				<a href="<?= "http://rfastofficial.github.io/?function=".$func; ?>"><?= $func; ?></a>

			</li>

	<?php

		endforeach; 

		echo "</ul>";

	?>

</div>



<?php 



include("./_footer.php"); 



?>