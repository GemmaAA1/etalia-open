define([
    'app',
    'text!app/templates/thread/comment/thumb.hbs',
    'app/model/thread/comment',
    'app/view/user/thumb',
    'app/view/thread/comment/form'
], function (App, template) {

    App.View.Thread.CommentThumb = App.Backbone.View.extend({

        tagName: 'div',
        className: 'thread-comment',

        template: App.Handlebars.compile(template),

        events: {
            "click .comment-edit": "onEditClick",
            "click .comment-delete": "onDeleteClick"
        },

        initialize: function () {
            this.listenTo(this.model, "sync", this.render);
        },

        onEditClick: function(e) {
            e.preventDefault();

            var $contents = this.$('.comment-content, .comment-info').hide();

            var form = App.View.Thread.CommentForm.create({
                model: this.model
            },{
                $target: this.$('[data-comment-edit-form]')
            });

            this.listenToOnce(form, 'validation_success', function() {
                form.model.save(null, {wait: true});
            });
            this.listenToOnce(form, 'cancel', function() {
                form.$el.replaceWith('<div data-comment-edit-form></div>');
                form.remove();
                $contents.show();
            });
        },

        onDeleteClick: function(e) {
            e.preventDefault();

            if (confirm('Are you sure you want to delete your comment ?')) {
                var that = this;
                this.model.destroy({
                    success: function() {
                        that.remove();
                    }
                });
            }
        },

        render: function () {
            App.log('CommentThumbView::render', this.model.get('id'));

            var attributes = App._.extend({}, this.model.attributes, {
                is_owner: this.model.isOwner(App.getCurrentUser())
            });

            this.$el.html(this.template(attributes));

            this.pushSubView(
                App.View.User.Thumb.create({
                    model: this.model.get('user')
                }, {
                    $target: this.$('[data-user-thumb]')
                })
            );

            return this;
        }
    });

    App.View.Thread.CommentThumb.create = function(options, createOptions) {
        options = options || {};
        if (!options.model) {
            throw 'options.model is expected to be an instance of App.Model.Comment';
        }
        if (!options.id) {
            options.id = 'thread-comment-' + options.model.get('id');
        }

        var thumb = new App.View.Thread.CommentThumb(options);

        if (createOptions) {
            App.View.create(thumb, createOptions);
        }

        return thumb;
    };

    return App.View.Thread.CommentThumb;
});
