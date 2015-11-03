$(document).ready(function() {

    $('.question').on('click', function () {
        $(this).next('.response').slideToggle(250);
    })

});
