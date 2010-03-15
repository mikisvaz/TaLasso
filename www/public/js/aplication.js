jQuery(document).ready(function(){ 

	//close all the content divs on page load 
	$('.mover').hide(); 
	
	// toggle slide 
	$('#slideToggle').click(function(){ 
		// by calling sibling, we can use same div for all demos 
		$(this).siblings('.mover').slideToggle(); 
	}); 
}); // --></script>
