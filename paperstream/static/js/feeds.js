$(document).ready(function() {

    toggleStamps();
    extendPaper();

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
    dislike();

    // Send Like/dislike ajax call
    like();
});


function toggleStamps(){
    $('.paper-list').on({
        mouseenter: function () {
            $(this).find('.stamps').show();
            $(this).find('.stamps-2').hide();
        },
        mouseleave: function () {
            $(this).find('.stamps').hide();
            $(this).find('.stamps-2').show();
        }
    });
}

function extendPaper(){
    $('.paper-list').click( function () {
        $(this).children('.compact').toggle();
        $(this).children('.extended').toggle();
        if ($(this).children('.extended').is(":visible")) {
            $(this).unbind('mouseenter mouseleave');
            $(this).children('.stamps').show();
        } else {
            $(this).on({
                mouseenter: function () {
                    $(this).find('.stamps').show();
                    $(this).find('.stamps-2').hide();
                },
                mouseleave: function () {
                    $(this).find('.stamps').hide();
                    $(this).find('.stamps-2').show();
                }
            });
        }
    }).children('.stamps').click(function(e) {
        e.stopPropagation();
    });
}

function dislike () {
    $('.dislike').on('click', function (event) {
        var dislike = this;
        var id = $(this).parents('span').attr('id');
        var url = $(location).attr('protocol') + '//' +
            $(location).attr('host') + '/user/paper/dislike';
        console.log(url);
        $.ajax({
            type: "POST",
            url: url,
            data: {pk: id},
            success: function (json) {
                var like = $(dislike).siblings(".like");
                var like2 = $(dislike).parent().siblings('.compact').find('.like');
                var dislike2 = $(dislike).parent().siblings('.compact').find('.dislike');
                console.log(like);
                console.log(like);
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
}

function like () {
    $('.like').on('click', function (event) {
        var like = this;
        var id = $(this).parents('span').attr('id');
        var url = $(location).attr('protocol') + '//' +
            $(location).attr('host') + '/user/paper/like';
        $.ajax({
            type: "POST",
            url: url,
            data: {pk: id},
            success: function (json) {
                var dislike = $(like).siblings(".dislike");
                var like2 = $(dislike).parent().siblings('.compact').find('.like');
                var dislike2 = $(dislike).parent().siblings('.compact').find('.dislike');
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
}