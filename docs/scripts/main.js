
function check_mail_ok (form,url,reset_pressed) {
	if(reset_pressed==true){
		form.reset();
		return false;
	}
	var name=form.elements["name"].value;
	var subject=form.elements["subject"].value;
	var receiver=form.elements["receiver"].value;
	var message=form.elements["message"].value;
	if(name=="" || subject=="selected" || receiver=="selected" || message==""){
		alert("All fields are mandatory. :)");
		return false;
	}
	form.action = url;
	return true;
}