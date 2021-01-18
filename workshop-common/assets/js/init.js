(function($){
  $(function(){

    $('.button-collapse').sideNav();
    $('.collapsible').collapsible();
    $('.parallax').parallax();
    $(".dropdown-button").dropdown({ 
            hover: true,
            constrainWidth: false,
            belowOrigin: true, 
    });
            

    $('select').material_select();
  }); // end of document ready
})(jQuery); // end of jQuery name space
