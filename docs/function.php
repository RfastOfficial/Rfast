<?php

include("./_head.php"); ?>

<div id='content'>

	<?php 
		$function_s = $sanitizer->text($input->get->function);
		$path=getcwd().'/functions/'.$function_s;
		$author=file_get_contents($path."/author");
		$description=file_get_contents($path."/description");
		echo "<b>Function name:</b> &nbsp;".$function_s."</br></br>";
		echo "<b>Description:</b> <p>&nbsp;&nbsp;".$description."<p></br>";
		echo "<b>Authors:</b> &nbsp;".$author."</br>";
	?>
</div>

<?php 

include("./_foot.php"); ?>