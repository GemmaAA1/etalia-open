
$(document).ready(function() {
    // Send Like/tick ajax call
    $('.like').on('click', like);

    // Send add-to-lib ajax call
    $('.add-to-library').on('click', add_to_lib);

    // Send trash paper from user library ajax call
    $('.trash').on('click', trash_paper);

});

function like (event) {
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
}

function add_to_lib (event) {
    var $add = $(this);
    $add.removeClass('add-to-library')
        .addClass('loading');
    var id = $(this).parents('ul').attr('id');
    var url = '/user/paper/add';
    $.ajax({
        type: 'POST',
        url: url,
        data: {pk: id},
        success: function (json) {
            $.each(json, function (key, value) {
                if (key == 'success') {
                    if (value) {  // success
                        console.log('add success');
                        $add.removeClass('loading')
                            .addClass('trash')
                            .parents('ul').find('.like').addClass('active');
                        $add.on('click', trash_paper)
                            .off('click', add_to_lib);
                    } else {  // fail
                        console.log('add failed');
                        $add.removeClass('loading')
                            .addClass('add-to-library');
                    }
                } else if (key == 'message') {
                    if (value) {
                        console.log(value);
                    }
                }
            });
        }
    });
    event.stopPropagation();
}

function trash_paper (event) {
    var $add = $(this);
    $add.removeClass('trash')
        .addClass('loading');
    var id = $(this).parents('ul').attr('id');
    var url = '/user/paper/add';
    $.ajax({
        type: 'POST',
        url: url,
        data: {pk: id},
        success: function (json) {
            $.each(json, function (key, value) {
                if (key == 'success') {
                    if (value) {
                        $add.removeClass('loading')
                            .addClass('add-to-library');
                        $add.on('click', add_to_lib)
                            .off('click', trash_paper);
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
}