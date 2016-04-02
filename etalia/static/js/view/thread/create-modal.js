define([
    'jquery',
    'backbone',
    'model/thread/thread',
    'backbone-forms',
    'backbone-modal'
], function ($, Backbone, ThreadModel) {

    var submitForm = function(e) {
        e.preventDefault();

        var errors = this.form.commit({ validate: true });
        if (errors) {
            console.log('Errors', errors);
            return;
        }

        this.form.model.on('sync', function(e) {
            console.log(e);
            console.log('Created !');
            //console.log(this.form.model.toJSON());

            // TODO redirect
        });

        console.log('Saving ...', this.form.model);
        this.form.model.save();
    };

    var ThreadCreateModal = Backbone.ModalView.extend({
        title: "Start a new thread",

        form: null,

        buttons: [{
            className: "btn-primary thread-create-modal-post",
            label: "Post"
        }, {
            className: "btn-default thread-create-modal-cancel",
            label: "Cancel",
            close: true
        }],

        events: {
            "click .thread-create-modal-post": "onPost",
            "click .thread-create-modal-cancel": "onCancel",
            "hidden.bs.modal": "onHidden"
        },

        postRender: function() {
            var model = ThreadModel.createNew();

            this.form = new Backbone.Form({
                model: model
            }).render();

            this.form.on('submit', submitForm);

            this.$body.append(this.form.el);

            return this;
        },

        onPost: submitForm,

        onCancel: function(e) {
            console.log("Cancel clicked");
        },

        onHidden: function(e) {
            console.log("Modal hidden");
        }
    });

    return ThreadCreateModal;
});
