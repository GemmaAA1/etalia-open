$(document).ready(function() {

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
        var element = this;
        var id = $(this).parents('li').attr('id');
        var url = $(location).attr('href') + 'dislike';
        $.ajax({
            type: "POST",
            url: url,
            data: {pk: id},
            success: function (json) {
                var button_dislike = $(element).children(":first");
                console.log($(button_dislike).attr('class'));
                var button_like = $(element).next().children(":first");
                console.log($(button_like).attr('class'));
                $.each(json, function (key, value) {
                    if (key == 'is_liked') {
                        if (value) {
                            $(button_like).removeClass('btn-default');
                            $(button_like).addClass('btn-info');
                        } else {
                            $(button_like).addClass('btn-default');
                            $(button_like).removeClass('btn-info');
                        }
                    } else if (key == 'is_disliked') {
                        if (value) {
                            $(button_dislike).removeClass('btn-default');
                            $(button_dislike).addClass('btn-info');
                        } else {
                            $(button_dislike).addClass('btn-default');
                            $(button_dislike).removeClass('btn-info');
                        }
                    }
                });
            }
        });
        event.preventDefault();
    });

    // Send Like/dislike ajax call
    $('.like').on('click', function (event) {
        var element = this;
        var id = $(this).parents('li').attr('id');
        var url = $(location).attr('href') + 'like';
        $.ajax({
            type: "POST",
            url: url,
            data: {pk: id},
            success: function (json) {
                var button_like = $(element).children(":first");
                console.log($(button_dislike).attr('class'));
                var button_dislike = $(element).prev().children(":first");
                console.log($(button_like).attr('class'));
                $.each(json, function (key, value) {
                    if (key == 'is_liked') {
                        if (value) {
                            $(button_like).removeClass('btn-default');
                            $(button_like).addClass('btn-info');
                        } else {
                            $(button_like).addClass('btn-default');
                            $(button_like).removeClass('btn-info');
                        }
                    } else if (key == 'is_disliked') {
                        if (value) {
                            $(button_dislike).removeClass('btn-default');
                            $(button_dislike).addClass('btn-info');
                        } else {
                            $(button_dislike).addClass('btn-default');
                            $(button_dislike).removeClass('btn-info');
                        }
                    }
                });
            }
        });
        event.preventDefault();
    });

    function update_like_dislike(key, value) {
        var $button_dislike = $(this).children('button');
        var $button_like = $(this).closest('span').children('button');
        if (key == 'is_liked') {
            if (value) {
                $button_like.removeClass('btn-default');
                $button_like.addClass('btn-info');
            } else {
                $button_like.addClass('btn-default');
                $button_like.removeClass('btn-info');
            }
        } else if (key == 'is_disliked') {
            if (value) {
                $button_dislike.removeClass('btn-default');
                $button_dislike.addClass('btn-info');
            } else {
                $button_dislike.addClass('btn-default');
                $button_dislike.removeClass('btn-info');
            }
        }
    }
});