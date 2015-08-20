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

    // To update seed-paper from modify-seed
    $('tr[seed-paper]').on('click', function (event) {
        var $tr = $(this);
        //console.log('form submitted!');
        $.ajax({
            type: "POST",
            url: $(location).attr('href'),
            data: {pk: $tr.attr('id')},
            success: function (json) {
                $.each(json, function (key, value) {
                    var $tr2 = $('#' + key);
                    $tr2.removeClass();
                    $tr2.addClass(value);
                });
            }
        });
        event.preventDefault();
    });

});