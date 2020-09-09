</div><!--/#main-->
	<!-- footer -->
	<footer id="footer">
		<div class="footer-content">
			<p>
			Powered by Rfast <!--&nbsp; / &nbsp; 
			?php 
			if($user->isLoggedin()) {
				// if user is logged in, show a logout link
				echo "<a href='{$config->urls->admin}login/logout/'>Logout ($user->name)</a>";
			} else {
				// if user not logged in, show a login link
				echo "<a href='{$config->urls->admin}'>Admin Login</a>";
			}
			?-->
			</p>
		</div>
	</footer>


</body></html>