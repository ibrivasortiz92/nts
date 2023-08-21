// DISPLAY SECTIONS
document.addEventListener("DOMContentLoaded", function () {

    // SECTION LINKS
    document.getElementById("setup-link").onclick = () => display_section("setup");
    document.getElementById("solution-link").onclick = () => display_section("solution");
    document.getElementById("results-link").onclick = () => display_section("results");

    // INITIAL PAGE
    display_section("setup");
  
});  
  
// DISPLAY SECTION FUNCTION
function display_section(section){

    // HIDE ALL SECTIONS
    document.getElementById("profile-section").style.display = "none";
    document.getElementById("main-content__layout").style.display = "none";
    document.getElementById("setup-section").style.display = "none";
    document.getElementById("solution-section").style.display = "none";
    document.getElementById("results-section").style.display = "none";

    // SHOW REQUESTED SECTION
    if(section == "profile") {
        document.getElementById("profile-section").style.display = "block";
    }
    else if (section == "setup") {
        document.getElementById("main-content__layout").style.display = "grid";
        document.getElementById("setup-section").style.display = "flex";
        setup_section_init();
    }
    else if (section == "solution") {
        document.getElementById("main-content__layout").style.display = "grid";
        document.getElementById("solution-section").style.display = "flex";
        solution_section_init();
    }
    else if (section == "results") {
        document.getElementById("main-content__layout").style.display = "grid";
        document.getElementById("results-section").style.display = "flex";
        results_section_init();
    }

}


// TOGGLE MENU
document.querySelector("#menu-toggle").onclick = function () {
    if (window.innerWidth > 800){
      	document.querySelector(".sidebar").classList.toggle("active");
    }
}

window.onresize = function() {
    if (document.getElementById("canvas_geom")) document.getElementById("canvas_geom").remove();
    if (document.getElementById("canvas_quad")) document.getElementById("canvas_quad").remove();
    if (document.getElementById("results__picture")) document.getElementById("results__picture").remove();
    document.getElementById("initial_img").style.display = "block";
}