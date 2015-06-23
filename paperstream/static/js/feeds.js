$(document).ready(function() {
    $(".list-group-item").each(function() {
        $(this).click(function (event) {
            var ufp_id;
            ufp_id = $(this).attr("id");
            if(!$(event.target).closest('#likes'+ufp_id).length &&
               !$(event.target).closest('#dislikes'+ufp_id).length){
                $("#abstract"+ufp_id).toggleClass('active');
            }
        });
    });
});