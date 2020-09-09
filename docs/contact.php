<?php

include("./_header.php"); ?>
<script src="scripts/main.js" type="text/javascript">var reset_pressed=false;</script>
<div id='content'>
	<?php
		if(isset($_POST['send'])){
			$name="Name: ".$_POST["name"])."\n\n";
			$email=$_POST["receiver"]);
			$message=$_POST["message"]);
			$subject=$_POST["subject"]);
			$ok=mail($email, $subject, $name.$message, "From: rfast");
			$url_to_img=$config->urls->templates."images/emoji.jpg";
			echo $ok ? "<b>The message was sent succesfully. </b><img src='$url_to_img' alt='Smiley' width='30' height='30'>" : 
					"<b>There was a problem, possibly with the network...try again... :(</b>";
		}
		else if(isset($_GET["emails"])){
			echo "<table class='table'>
					  <thead>
					    <tr>
					      <th></th>
					      <th>Name</th>
					      <th>Email</th>
					    </tr>
					  </thead>
					  <tbody>
					    <tr>
					      <th scope='row'>1</th>
					      <td>Manos Papadakis</td>
					      <td>papadakm95@gmail.com</td>
					    </tr>
					    <tr>
					      <th scope='row'>2</th>
					      <td>Michail Tsagris</td>
					      <td>mtsagris@yahoo.gr</td>
					    </tr>
					  </tbody>
					</table>";
		}else if(isset($_GET["message"])){
			echo "<form onsubmit=\"return check_mail_ok(this,$page->url,reset_pressed);\" class='contact' method='post'>
				    <fieldset class='contact-inner'>
				      <p class='contact-input'>
				        <input type='text' name='name' placeholder='Your Name...' autofocus>
				      </p>

					  <p class='contact-input'>
				        <input type='email' name='email' placeholder='Your email...' autofocus>
				      </p>

				      <p class='contact-input'>
				        <label for='select' class='select'>
				          <select name='subject' id='select'>
				            <option value='selected'>Choose Subject...</option>
				            <option value='I have a suggestion'>I have a suggestion</option>
				            <option value='I found a bug'>I found a bug</option>
				            <option value='Other'>Other</option>
				          </select>
				        </label>
				      </p>

				      

				      <p class='contact-input'>
				        <textarea name='message' placeholder='Your Message...'></textarea>
				      </p>

				      <p class='contact-submit'>
				        <input type='submit' name='send' onclick='reset_pressed=false;' value='Send Message'>
				      </p>
				      <p class='contact-reset'>
				        <input type='submit' onclick='reset_pressed=true;' name='reset_b' value='reset'>
				      </p>
				    </fieldset>
				</form>";
		}
		?>

</div>
<?php 

include("./footer.php"); ?>