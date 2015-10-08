$(document).ready(function() {

    $( ".paper-list" ).hover(function() {
        var $test = $(this).find('.stamps');
        $(this).find('.stamps').toggle();
        $(this).find('.stamps-2').toggle();
    });

    $(".list-group-item").each(function() {
        $(this).click(function (event) {
            var ufp_id;
            ufp_id = $(this).attr("id");
            if(!$(event.target).closest('.like').length &&
               !$(event.target).closest('.dislike').length &&
               !$(event.target).closest('#paper-title'+ufp_id).length){
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

    // Send Like/dislike ajax call
    $('.dislike').on('click', function (event) {
        var dislike = this;
        var id = $(this).parents('span').attr('id');
        var url = $(location).attr('href') + 'dislike';
        $.ajax({
            type: "POST",
            url: url,
            data: {pk: id},
            success: function (json) {
                var like = $(dislike).siblings(".like");
                var like2 = $(dislike).parents(":first").siblings(".stamps-2").children('.like');
                var dislike2 = $(dislike).parents(":first").siblings(".stamps-2").children('.dislike');
                $.each(json, function (key, value) {
                    if (key == 'is_liked') {
                        if (value) {
                            $(like).addClass('active');
                            $(like2).addClass('active');
                        } else {
                            $(like).removeClass('active');
                            $(like2).removeClass('active');
                        }
                    } else if (key == 'is_disliked') {
                        if (value) {
                            $(dislike).addClass('active');
                            $(dislike2).addClass('active');
                        } else {
                            $(dislike).removeClass('active');
                            $(dislike2).removeClass('active');
                        }
                    }
                });
            }
        });
        event.preventDefault();
    });

    // Send Like/dislike ajax call
    $('.like').on('click', function (event) {
        var like = this;
        var id = $(this).parents('span').attr('id');
        var url = $(location).attr('href') + 'like';
        $.ajax({
            type: "POST",
            url: url,
            data: {pk: id},
            success: function (json) {
                var dislike = $(like).siblings(".dislike");
                var like2 = $(dislike).parents(":first").siblings(".stamps-2").children('.like');
                var dislike2 = $(dislike).parents(":first").siblings(".stamps-2").children('.dislike');
                console.log(like2);
                console.log(dislike2);
                $.each(json, function (key, value) {
                    if (key == 'is_liked') {
                        if (value) {
                            $(like).addClass('active');
                            $(like2).addClass('active');
                        } else {
                            $(like).removeClass('active');
                            $(like2).removeClass('active');
                        }
                    } else if (key == 'is_disliked') {
                        if (value) {
                            $(dislike).addClass('active');
                            $(dislike2).addClass('active');
                        } else {
                            $(dislike).removeClass('active');
                            $(dislike2).removeClass('active');
                        }
                    }
                });
            }
        });
        event.preventDefault();
    });

});