define([
    'app',
    'app/view/modal',
    'app/view/thread/list',
    'app/view/thread/form-create'
], function (App) {

    var list = new App.View.Thread.List({
            el: '#thread-list-container'
        });

    list.render();

});
