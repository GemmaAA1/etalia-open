$(document).ready(function() {

    toggleStamps();
    extendPaper();

    $(".list-group-item").each(function() {
        $(this).click(function (event) {
            var ufp_id;
            ufp_id = $(this).attr("id");
            if(!$(event.target).closest('.like').length &&
               !$(event.target).closest('.tick').length &&
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

    // Send Like/tick ajax call
    tick();

    // Send Like/tick ajax call
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
        //$(this).children('.extended').slideToggle();
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

function tick () {
    $('.tick').on('click', function (event) {
        var tick = this;
        var id = $(this).parents('span').attr('id');
        var url = $(location).attr('protocol') + '//' +
            $(location).attr('host') + '/user/paper/tick';
        console.log(url);
        $.ajax({
            type: "POST",
            url: url,
            data: {pk: id},
            success: function (json) {
                var like = $(tick).siblings(".like");
                var like2 = $(tick).parent().siblings('.compact').find('.like');
                var tick2 = $(tick).parent().siblings('.compact').find('.tick');
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
                    } else if (key == 'is_ticked') {
                        $(tick).parents('.paper-list').slideUp(250);
                        if (value) {
                            $(tick).addClass('active');
                            $(tick2).addClass('active');
                        } else {
                            $(tick).removeClass('active');
                            $(tick2).removeClass('active');
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
                var tick = $(like).siblings(".tick");
                var like2 = $(tick).parent().siblings('.compact').find('.like');
                var tick2 = $(tick).parent().siblings('.compact').find('.tick');
                $.each(json, function (key, value) {
                    if (key == 'is_liked') {
                        if (value) {
                            $(like).addClass('active');
                            $(like2).addClass('active');
                        } else {
                            $(like).removeClass('active');
                            $(like2).removeClass('active');
                        }
                    } else if (key == 'is_ticked') {
                        if (value) {
                            $(tick).addClass('active');
                            $(tick2).addClass('active');
                        } else {
                            $(tick).removeClass('active');
                            $(tick2).removeClass('active');
                        }
                    }
                });
            }
        });
        event.preventDefault();
    });
}