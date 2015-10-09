
$(document).ready(function() {

    // Send Like/dislike ajax call
    $('.dislike').on('click', function (event) {
        var dislike = this;
        var id = $(this).parent('div').attr('id');
        var url = $(location).attr('protocol') + '//' + $(location).attr('host') + '/user/paper/dislike';
        $.ajax({
            type: "POST",
            url: url,
            data: {pk: id},
            success: function (json) {
                var like = $(dislike).siblings(".like");
                $.each(json, function (key, value) {
                    if (key == 'is_liked') {
                        if (value) {
                            $(like).addClass('active');
                        } else {
                            $(like).removeClass('active');
                        }
                    } else if (key == 'is_disliked') {
                        if (value) {
                            $(dislike).addClass('active');
                        } else {
                            $(dislike).removeClass('active');
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
        var id = $(this).parent('div').attr('id');
        var url = $(location).attr('protocol') + '//' + $(location).attr('host') + '/user/paper/like';
        $.ajax({
            type: "POST",
            url: url,
            data: {pk: id},
            success: function (json) {
                var dislike = $(like).siblings(".dislike");
                $.each(json, function (key, value) {
                    if (key == 'is_liked') {
                        if (value) {
                            $(like).addClass('active');
                        } else {
                            $(like).removeClass('active');
                        }
                    } else if (key == 'is_disliked') {
                        if (value) {
                            $(dislike).addClass('active');
                        } else {
                            $(dislike).removeClass('active');
                        }
                    }
                });
            }
        });
        event.preventDefault();
    });
});
