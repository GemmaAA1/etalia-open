
$(document).ready(function() {

    // Send Like/tick ajax call
    like();

    //
    add_to_lib();
});


function like () {
    $('.like').on('click', function (event) {
        var like = this;
        var id = $(this).parents('ul').attr('id');
        var url = $(location).attr('protocol') + '//' +
            $(location).attr('host') + '/user/paper/like';
        $.ajax({
            type: 'POST',
            url: url,
            data: {pk: id},
            success: function (json) {
                $.each(json, function (key, value) {
                    if (key == 'is_liked') {
                        if (value) {
                            $(like).addClass('active');
                        } else {
                            $(like).removeClass('active');
                        }
                    }
                });
            }
        });
        event.preventDefault();
    });
}

function add_to_lib () {
    $('.add-to-library').on('click', function (event) {
        var add = this;
        var id = $(this).parents('ul').attr('id');
        var url = $(location).attr('protocol') + '//' +
            $(location).attr('host') + '/user/paper/add';
        $.ajax({
            type: 'POST',
            url: url,
            data: {pk: id},
            success: function (json) {
                $.each(json, function (key, value) {
                    if (key == 'success') {
                        if (value) {
                            $(add).parents('ul').find('.like').addClass('active');
                            $(add).replaceWith('<span class="add-to-library active"></a>');
                        } else {}
                    } else if (key == 'message') {
                        if (value) {
                            console.log(value);
                        }
                    }
                });
            }
        });
        event.stopPropagation();
    });
}