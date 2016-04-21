define([
    'jquery',
    'app',
    'app/view/modal',
    'app/view/thread/list',
    'app/view/thread/form-create',
    'bootstrap'
], function ($, App) {

    // Thread list
    var list = new App.View.Thread.List({
            el: '#thread-list-container'
        });

    list.render();


    var form, modal;
    $('#thread-create-modal').on('click', function () {
        form = App.View.Thread.CreateForm.create();

        modal = new App.View.Modal({
            title: 'Start a new thread',
            content: form,
            footer: false
        });

        form.once('validation_success', function () {
            form.model.save(null, {
                success: function () {
                    list.collection.add(form.model);
                    modal.close();
                    list.openDetail(form.model);
                },
                error: function () {
                    // TODO
                }
            });
        });

        form.once('cancel', function () {
            modal.close();
        });

        modal.once('hidden', function () {
            form = null;
            modal = null;
        });

        modal.render();
    });

});
