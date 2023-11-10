function header() {
  
  // add the header structure as the firstChild of div.reveal-header
  function add_header() {
    let header = document.querySelector("div.reveal-header");
    let reveal = document.querySelector(".reveal");
    reveal.insertBefore(header, reveal.firstChild);
    
    logo_img_left = document.querySelector('.header-logo-left img');
    if (logo_img_left.getAttribute('src') == null) {
      if (logo_img_left?.getAttribute('data-src') != null) {
        logo_img_left.src = logo_img_left?.getAttribute('data-src') || "";
        logo_img_left.removeAttribute('data-src'); 
      };
    };
    
    logo_img_right = document.querySelector('.header-logo-right img');
    if (logo_img_right.getAttribute('src') == null) {
      if (logo_img_right?.getAttribute('data-src') != null) {
        logo_img_right.src = logo_img_right?.getAttribute('data-src') || "";
        logo_img_right.removeAttribute('data-src'); 
      };
    };
  };
  
  // dynamically changing the header
  function change_header(dheader, cheader, ctext) {
    // dhead => dynamic header
    // chead => constant header
    // ctext => contstant header_text inner html
    if (dheader !== null) {
      cheader.innerHTML = dheader.innerHTML;  
    } else {
      cheader.innerHTML = ctext;
    };
  };
  
  if (Reveal.isReady()) {
    add_header();
  }; 
};

window.addEventListener("load", (event) => {
  header();
});
