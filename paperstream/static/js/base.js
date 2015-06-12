$(document).ready(function() {
    $(".likes").each( function() {
        $(this).click( function(event) {
            var fnp_id;
            fnp_id = $(this).attr("data-fnpid");
            $.getJSON('/feed/likes', {feednewpaper_id: fnp_id},
                function(json) {

                var fnp_id = json.fnp_id;

                if (json.is_liked) {
                    $("#" + fnp_id).fadeTo(200, 1);
                    $("#dislikes" + fnp_id).children().removeClass("btn-info");
                    $("#dislikes" + fnp_id).children().addClass("btn-default");
                    $("#likes" + fnp_id).children().addClass("btn-info");
                    $("#likes" + fnp_id).children().removeClass("btn-default");
                }
                else {
                    $("#likes" + fnp_id).children().addClass("btn-default");
                    $("#likes" + fnp_id).children().removeClass("btn-info");
                }
            });
        });
    });
});
