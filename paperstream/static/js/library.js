
$(document).ready(function() {

    // Send Like/tick ajax call
    like();
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